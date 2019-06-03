       function ampNonresonantHeavyImC4MM()
           implicit none
           type(dd_complex) :: ampNonresonantHeavyImC4MM

           type(dd_complex) :: result

      result = ddreal(0)

           ampNonresonantHeavyImC4MM = result
       end function ampNonresonantHeavyImC4MM

       function ampNonresonantHeavyImC4MP()
           implicit none
           type(dd_complex) :: ampNonresonantHeavyImC4MP
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           type(dd_complex) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           type(dd_complex) ::  t31,t32,t33,t34,t35,t36,t37,t38,t4,t5,t6,t7
           type(dd_complex) ::  t8,t9

           type(dd_complex) :: result

      t1 = intHLs250000x0121D2eps0()
      t2 = intHLs250000x0121D2eps1() * epinv
      t3 = (t2 + t1) * s25 + epinv2
      t4 = (gb)**2
      t5 = propZ25 * s25
      t6 = t5 * (gw)**2
      t7 = s25 + s26 + s56
      t8 = s16 + s26 - s34 + s56
      t9 = ddreal(2) * propW16
      t10 = ddreal(6) * propW16
      t11 = ddreal(3) * t6
      t12 = intHLs250000x0112D2eps1()
      t13 = intHLs250000x0113D4eps1()
      t14 = intHLs250000x0122D4eps1()
      t15 = ddreal(2) * t13
      t16 = (-t15 + t2 + t12 - t14 + t1) * s25 + epinv2
      t17 = t9 * t8 - ddreal(1)
      t18 = s25 + s16 - s34
      t19 = s26 + s56
      t20 = t18 + ddreal(2) * t19
      t21 = (s34)**2
      t22 = (s16)**2
      t23 = (s25)**2 + t21 + t22
      t24 = (s56)**2
      t21 = s25 * (s16 * s25 - s25 * s34 + t21 + t22)
      t22 = s12 * t7
      t25 = s26 * propW16
      t26 = t25 * ((s25 + s56) * s16 + (-s25 - s56) * s34 + s56 * (s25 +
     & s26))
      t27 = (s26)**2
      t28 = s16 * s34
      t28 = (((-s16 + s34) * s56 + t28) * s25 - s26 * t27 + (t28 - t24) 
     &* s26 + t28 * s56) * propW16 - s45 * t19 - s46 * t19
      t29 = t27 * propW16 * (s25 + s16 - s34)
      t30 = t10 * t8 + t5 - ddreal(4)
      t31 = -ddreal(4) * s25 + ddreal(2) * s34 - ddreal(6) * t19 - s16 *
     & (t5 + ddreal(2)) + t5 * (s25 + s34)
      t32 = -ddreal(1) + t5
      t33 = t32 * s26
      t32 = t32 * s56
      t34 = s25 * (-ddreal(1) + (s25 + s45 + s46) * propZ25)
      t35 = -s45 - s46
      t36 = t4 * t30
      t37 = t11 * t17 + t36
      t38 = (t2 + t1 - t14) * t4
      t9 = (t4 * (t10 * t7 + t5 + ddreal(2)) + t11 * (t9 * t7 + ddreal(1
     &))) * t8
      t1 = t16 * (-struc2MP * (t11 * t20 - t4 * t31) + t7 * t37 * (struc
     &3MP + struc8MP)) - t7 * (epinv2 * (t4 * (-t5 + ddreal(4)) + t11) +
     & s25 * (ddreal(6) * t6 * t13 * t17 + ddreal(6) * (t4 * t8 + t6 * t
     &8) * t14 * propW16 - t12 * t37 - t38 * t5 + ddreal(3) * t6 * (t2 +
     & t1 - t14) + ddreal(4) * t38 + t36 * t15)) * struc5MP + t9 * (-t15
     & + t12 - t14) * s25 * struc7MP
      t2 = ddreal(1) / t7
      t10 = ddreal(1) / t8
      result = -ddreal(2) / ddreal(3) * im * propW34 * (t9 * t3 * (struc
     &15MP + struc41MP + struc57MP) - t16 * (t4 * (ddreal(12) * t26 + dd
     &real(9) * t29 - ddreal(6) * t28 + ddreal(3) * (propW16 * t23 + s25
     &) * s56 + ddreal(3) * t24 * (propW16 * t18 + ddreal(1)) + ddreal(3
     &) * t25 * t23 + ddreal(3) * t21 * propW16 + s13 * t31 + (t34 + t33
     & + t32) * s16 + (s26 * (propZ25 * t19 - ddreal(1)) + t5 * (s26 - s
     &45 - s46)) * s25 + (-t34 - t33 - t32) * s34 - t27 + t22 * t30 - dd
     &real(2) * s16 * t35 + ddreal(2) * s26 * s56 + ddreal(2) * s34 * t3
     &5 - ddreal(4) * s25 * t35) + t11 * (ddreal(4) * t26 + (s26 * t23 +
     & s56 * t23 + t18 * t24 + t21) * propW16 - s13 * t20 + s45 * t18 + 
     &s46 * t18 + s56 * (s25 + s26) + t24 + t22 * t17 - ddreal(2) * t28 
     &+ ddreal(3) * t29)) * struc1MP - ddreal(2) * t1 + ddreal(6) * t8 *
     & t7 * (t6 + t4) * t3 * propW16 * (struc14MP + struc40MP + struc56M
     &P)) / (ecossin)**2 / s25 * t2 * t10

           ampNonresonantHeavyImC4MP = result
       end function ampNonresonantHeavyImC4MP

       function ampNonresonantHeavyImC4PM()
           implicit none
           type(dd_complex) :: ampNonresonantHeavyImC4PM
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           type(dd_complex) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           type(dd_complex) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t5,t6
           type(dd_complex) ::  t7,t8,t9

           type(dd_complex) :: result

      t1 = s12 + s13
      t2 = ddreal(2) * s13
      t3 = s25 + s16 - s34
      t4 = ddreal(2) * s26
      t5 = t4 + t3
      t6 = s25 + s26
      t7 = s13 - s45
      t8 = s25 + s26 + s56
      t9 = s16 + s26 - s34 + s56
      t10 = s13 * s16
      t11 = propW16 * (t3 + ddreal(2) * s12 + ddreal(2) * s26) * t8 * t9
      t12 = t1 * s25
      t13 = t6 * s46
      t14 = intHLs250000x0121D2eps0()
      t15 = intHLs250000x0112D2eps1()
      t16 = s13 - s45 - s46
      t17 = -s45 - s46
      t18 = s16 - s34
      t19 = t18 * s13 + (t2 + s12) * s26 + (-s25 + s12 - s26 + ddreal(2)
     & * t16 - s56) * s56 + t5 * t17 - t11 + t12
      t20 = intHLs250000x0113D4eps1()
      t21 = intHLs250000x0122D4eps1()
      t1 = s16 + ddreal(4) * t1
      t22 = ddreal(6)
      t23 = t22 * s13
      t24 = ddreal(4) * s12
      t25 = ddreal(3)
      t26 = t25 * s26
      t27 = ddreal(2) * s25 + t18 + t26
      t28 = t25 * s25
      t29 = s25 - s16 + s34
      t30 = (s12 + s13 + s16) * s25
      t31 = (s56)**2
      t4 = ddreal(2) * t17 * t27 - t25 * (t11 + t31) + ddreal(2) * t10 +
     & (propZ25 * ((-s12 - s16 - s25 - s26) * s26 + (s25 - s13 + s26) * 
     &s34 - (s12 + s16 + s26 - s34) * s56 - t17 * t29 + t10 - t30) + t1)
     & * s25 + (t24 + t23 + s25 + s16 + s26) * s26 - (t2 + t6) * s34 + (
     &t16 * t22 + t18 + t24 - t28 - t4) * s56
      t6 = intHLs250000x0121D2eps1() * epinv
      t16 = (t6 + t14) * s25 + epinv2
      t17 = (gb)**2
      t32 = s25 * (t15 - t21)
      t33 = propZ25 * s25
      t34 = t33 * (gw)**2
      t35 = t33 + ddreal(2)
      t36 = t34 * t25
      t37 = s26 + s56
      t38 = ddreal(2) * propW16
      t39 = t22 * propW16
      t3 = ((-ddreal(2) * t20 + t6 + t14 + t15 - t21) * s25 + epinv2) * 
     &(-(t17 * (t22 * t37 + ddreal(4) * s25 - ddreal(2) * s34 + s16 * t3
     &5 - t33 * (s25 + s34)) + t36 * (t3 + ddreal(2) * t37)) * struc2PM 
     &+ t9 * (t17 * (t39 * t8 + t33 + ddreal(2)) + t36 * (t38 * t8 + ddr
     &eal(1))) * struc7PM + (t17 * (t39 * t9 + t33 - ddreal(4)) + t36 * 
     &(t38 * t9 - ddreal(1))) * t8 * (struc3PM + struc5PM + struc8PM))
      t15 = ddreal(1) / t9
      t8 = ddreal(1) / t8
      result = ddreal(2) / ddreal(3) * im * propW34 * ((t25 * t34 * (t16
     & * ((t2 + s16 + s12 + s26) * s26 - (s13 + s26) * s34 - t5 * s45 + 
     &(-s25 + s12 - s46 + ddreal(2) * t7 - s56) * s56 + t10 - t11 + t12 
     &- t13) + t32 * t19) + t17 * (t16 * (-ddreal(2) * s45 * t27 - t25 *
     & (-(s26)**2 + t11 + t31) + ddreal(2) * t10 + (propZ25 * (-(s12 + s
     &25) * s26 + (s25 - s13) * s34 + s45 * t29 - (s12 + s16 - s34 - s46
     &) * s56 + t10 + t13 - t30) + t1) * s25 + (s16 * t25 + s25 + t23 + 
     &t24) * s26 - (t26 + t2 + s25) * s34 + (t22 * t7 + ddreal(4) * s12 
     &- ddreal(4) * s46 + t18 - t28) * s56 - ddreal(4) * t13) + t32 * t4
     &) - t20 * s25 * (t34 * t22 * t19 + ddreal(2) * t17 * t4)) * struc1
     &PM + (t17 * t35 + t36) * ((t6 + t14) * s25 + epinv2) * t9 * (struc
     &27PM + struc43PM - struc57PM) - ddreal(2) * t3) / (ecossin)**2 / s
     &25 * t8 * t15

           ampNonresonantHeavyImC4PM = result
       end function ampNonresonantHeavyImC4PM

       function ampNonresonantHeavyImC4PP()
           implicit none
           type(dd_complex) :: ampNonresonantHeavyImC4PP

           type(dd_complex) :: result

      result = ddreal(0)

           ampNonresonantHeavyImC4PP = result
       end function ampNonresonantHeavyImC4PP

       function ampNonresonantHeavyImC7MM()
           implicit none
           type(dd_complex) :: ampNonresonantHeavyImC7MM

           type(dd_complex) :: result

      result = ddreal(0)

           ampNonresonantHeavyImC7MM = result
       end function ampNonresonantHeavyImC7MM

       function ampNonresonantHeavyImC7MP()
           implicit none
           type(dd_complex) :: ampNonresonantHeavyImC7MP
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           type(dd_complex) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           type(dd_complex) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           type(dd_complex) ::  t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52
           type(dd_complex) ::  t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63
           type(dd_complex) ::  t64,t65,t66,t67,t68,t69,t7,t70,t71,t72,t73,t74
           type(dd_complex) ::  t75,t76,t77,t78,t79,t8,t80,t81,t82,t83,t84,t85
           type(dd_complex) ::  t86,t87,t88,t89,t9,t90

           type(dd_complex) :: result

      t1 = intHLs250000x0112D2eps0()
      t2 = s25 + s16 - s12
      t3 = ddreal(2)
      t4 = t3 * s13
      t5 = s25 + s16 - s34
      t6 = t3 * s26
      t7 = s25 - s13 + s45 + s46
      t8 = s16 - s12 - s34
      t9 = ddreal(3)
      t10 = t9 * s26
      t11 = s12 + s26
      t12 = t3 * t11
      t13 = s25 + s26 + s56
      t14 = s16 + s26 - s34 + s56
      t15 = -s45 - s46
      t16 = s13 * s16
      t17 = (s12 + s13 - s16) * s25
      t18 = (-s13 + s25 + s26) * s34
      t19 = t15 * (t6 + t5)
      t20 = propW16 * t13
      t21 = t20 * (t12 + s25)
      t22 = t21 * t14
      t23 = (s56)**2
      t24 = t3 * t23 - (t4 - t2 - s26) * s26 + (t3 * t7 + t10 + t8) * s5
     &6 - t16 - t17 - t18 - t19 + t22
      t25 = s12 + s13
      t26 = t3 * t25
      t27 = t26 - s16
      t28 = s16 + s25
      t29 = t3 * s12
      t30 = t9 * s13
      t31 = t29 + t30 - t28
      t32 = s16 - s34
      t33 = t3 * s25
      t34 = s12 + s16 + s25
      t35 = -s16 + s25 + s34
      t36 = s12 + s16 + s26 - s34
      t37 = (s12 + s13 + s16) * s25
      t38 = t15 * t35
      t39 = (-t34 - s26) * s26 - s56 * t36 + t16 + t18 - t37 - t38
      t40 = (s26)**2
      t7 = (-t7 * t9 - ddreal(4) * s26 + t29 - t32) * s56
      t41 = t15 * (t33 + t10 + t32)
      t42 = propZ25 * t39
      t43 = ddreal(6) * t23
      t44 = t22 * t9
      t45 = (gb)**2
      t46 = (t45)**2
      t47 = (gw)**2
      t48 = (t47)**2
      t49 = t9 * s25
      t50 = t49 * t47 * propZ25
      t51 = intHLs250000x0113D4eps0()
      t52 = ddreal(4) * t25
      t53 = s25 + s26
      t54 = s13 - s45 - s46
      t55 = ddreal(9)
      t56 = ddreal(4) * s12
      t57 = ddreal(6) * t54
      t58 = t3 * t54
      t59 = t20 * (t3 * t36 + s25)
      t60 = t32 * s13
      t61 = ddreal(10)
      t62 = t49 * propZ25
      t61 = (-ddreal(8) * s12 - ddreal(12) * s13) * s26 + t14 * (ddreal(
     &6) * t20 * (t12 + t35) + t62 * (-t58 + t59)) + ddreal(18) * t23 - 
     &t3 * ((-ddreal(5) * s16 + t52) * s25 + (-s25 * t55 - ddreal(14) * 
     &s26 - ddreal(5) * t32 + t56 + t57) * s56) + t61 * ((t28 + s26) * s
     &26 - s34 * t53) - ddreal(4) * t41 - ddreal(4) * t60
      t63 = propZ25 * s25
      t64 = t63 * t45
      t65 = t64 * t39
      t66 = t63 * t55
      t67 = t66 * t48
      t68 = t67 * t24
      t69 = ddreal(6) * s13
      t59 = t59 * t14
      t70 = (s13 + s25 + s26) * s34
      t71 = (t4 + t34 + s26) * s26 + (t58 + t36) * s56 + t16 + t19 + t37
     & - t59 - t70
      t34 = t3 * (-(t3 * t53 + s13) * s34 + t16 + t41) - t9 * (t63 * t71
     & + t59) + ddreal(4) * t37 + ddreal(4) * t40 + (ddreal(4) * t34 + t
     &69) * s26 + (ddreal(4) * t36 + t57) * s56
      t36 = intHLs250000x0111D2eps0()
      t37 = ddreal(4) * t32
      t2 = -t3 * ((t26 + s16) * s25 + (t29 + t30 + t28 + s26) * s26 + (t
     &54 * t9 - s56 + t29 - t35) * s56 + t16 + t41 - t70) + t20 * (ddrea
     &l(6) * t11 + t49 + t37) * t14 + t63 * (t3 * (t20 * t32 * t14 - t23
     &) - (t2 + s26) * s26 + (-t33 - t10 - t8) * s56 - t16 + t17 + t38 +
     & t70)
      t8 = intHLs250000x0112D2eps1()
      t10 = s25 * t27 + (t31 - s26) * s26 + t16 + t18 + t41 + t7
      t11 = s26 + s56
      t17 = t3 * (s25 - s13 + s26 + s45 + s46 + s56) + t21
      t26 = s25 + s26 + s45 + s46 + s56
      t29 = ddreal(1) - epinv
      t30 = t63 * t48 * t24
      t38 = intHLs250000x0113D4eps1()
      t54 = -t3 + epinv
      t59 = intHLs250000x0122D4eps1()
      t70 = intHLs250000x0111D2eps1()
      t12 = t20 * (t12 + t5) * t14
      t72 = t9 * epinv
      t73 = t63 * t47
      t74 = t45 * intHLs250000x0122D4eps0()
      t2 = (-t36 * (t2 * t45 - t50 * t71) - t70 * (t45 * (epinv * t2 + t
     &3 * ((t42 + t52 + s16) * s25 + (t56 + t69 + t28 + s26) * s26 - (t4
     & + t53) * s34 + (-t6 - t49 + t56 + t57 + t32) * s56) - ddreal(6) *
     & t12 + ddreal(4) * t16 - ddreal(6) * t23 + ddreal(4) * t41) - t73 
     &* (-ddreal(6) * s25 * t25 - ddreal(6) * (t4 + s12) * s26 + ddreal(
     &6) * (-t58 + s25 - s12 + s26 + s56) * s56 - ddreal(6) * t19 - ddre
     &al(6) * t60 + ddreal(6) * t12 + t72 * t71)) + s25 * t1 * (-t45 * (
     &t3 * (t16 + t18 + t7 - t40 + t41) + t3 * t31 * s26 + (t3 * t27 + t
     &42) * s25 - t43 - t44) + t50 * t24)) * t47
      t4 = t3 * t32
      t6 = ddreal(6) * t11
      t7 = (propZ25 * t35 - ddreal(4)) * s25 - t4 - t6
      t5 = t3 * t11 + t5
      t12 = ddreal(1) / (gw)**2
      t16 = t45 * t12
      t18 = t5 * t9
      t19 = -t45 * t7 + t50 * t5
      t24 = propZ25 * t14
      t25 = t47 * t5
      t27 = t25 * t55
      t28 = t46 * t12
      t31 = t3 * ((t24 * t9 + ddreal(4)) * s25 + t4 + t6)
      t40 = t64 * t35
      t24 = (ddreal(6) * t24 + ddreal(8)) * s25
      t41 = ddreal(12) * t11
      t42 = t40 * t12
      t49 = ddreal(1) / s25
      t52 = t54 * t70
      t56 = t38 * t54
      t57 = t3 * propW16 * t14 - ddreal(1)
      t58 = ddreal(6) * propW16 * t14
      t60 = -ddreal(4) + t58 + t63
      t69 = t45 * t60
      t70 = im * propW34
      t71 = t9 * t47 * propZ25
      t75 = ddreal(1) - t63
      t76 = -t58 * t75 - t62 + ddreal(4)
      t77 = t63 + t3
      t78 = t47 * t57
      t79 = t3 * (t9 * propW16 * t77 * t14 - ddreal(4))
      t80 = t64 * t12
      t81 = ddreal(6) * t63
      t82 = t80 * t29
      t83 = ddreal(6) * t20
      t75 = t83 * t75
      t84 = t20 * t3 + ddreal(1)
      t83 = t83 + t63 + t3
      t85 = t45 * t83
      t86 = t20 * t77 + t63
      t87 = t47 * t84
      t88 = t3 * t45 * (t86 * t9 + t3)
      t89 = ddreal(1) / t14
      t90 = ddreal(1) / t13
      t4 = (-t3 * (-t51 * (-t31 * t45 + t63 * (t28 * t35 - t27)) + t56 *
     & (t67 * t5 - t45 * (-t31 * t47 + t40)) * t12) + t9 * (t1 * t19 + t
     &49 * (t45 * t36 * t7 - t52 * t19)) + t59 * (t25 * t66 - t45 * (t42
     & * t29 + epinv * ((t18 * propZ25 - ddreal(4)) * s25 - t4 - t6) - t
     &24 - t37 - t41)) + t8 * (-t55 * t73 * t5 * t29 + t45 * (-t72 * t7 
     &- t24 - t37 - t41 + t42)) + t74 * (s25 * ((t16 * t35 - t18) * prop
     &Z25 + ddreal(4)) + t4 + t6) - t27 * t36 * propZ25) * struc2MP * t8
     &9
      t5 = t70 * t90
      t2 = t5 * (t3 * s25 * (t51 * (t68 - t45 * (-t47 * t61 + t65)) + t3
     &8 * (t30 * (epinv * t55 - ddreal(18)) - t45 * (-t47 * (epinv * t61
     & - t14 * (ddreal(6) * t63 * t17 + ddreal(12) * t21) + ddreal(8) * 
     &t10 - ddreal(24) * t23) + t65 * t54))) - t9 * t2 - (-t8 * (t55 * t
     &30 * t29 - t45 * (-t47 * (-ddreal(8) * s12 * t13 - ddreal(8) * s25
     & * t15 + t9 * (epinv * (t10 * t3 + t63 * t39 - t43 - t44) + t63 * 
     &t14 * t17) - ddreal(12) * s26 * t15 + ddreal(12) * (s25 + s45 + s4
     &6) * s56 + ddreal(12) * t23 - ddreal(4) * s13 * (t11 * t9 + t32 + 
     &t33) + ddreal(4) * t26 * s16 + ddreal(4) * s26 * t53 - ddreal(4) *
     & t26 * s34 + ddreal(16) * s26 * s56 + ddreal(6) * t22) + t65)) - t
     &59 * (-t68 + t45 * (t47 * (epinv * t34 - t14 * (t62 * t17 + ddreal
     &(6) * t21) + ddreal(4) * t10 - ddreal(12) * t23) + t65 * t29)) + t
     &74 * (-t34 * t47 + t65)) * s25) * struc1MP * t12 * t49
      t2 = t89 * ((struc3MP + struc5MP + struc8MP) * (t3 * t70 * (t59 * 
     &(-t45 * (t82 - epinv * t76 + propW16 * t14 * (ddreal(12) + t81) - 
     &ddreal(8)) - t78 * t66) + t8 * (t55 * t73 * t57 * t29 + t45 * (t58
     & * t77 - t72 * t60 + t80 - ddreal(8))) + t74 * (t47 * t76 + t64) *
     & t12) + ddreal(6) * t70 * (t36 * (t69 * t49 + t71 * t57) + t52 * (
     &t50 * t57 + t69) * t49) + ddreal(4) * t70 * (t51 * (t79 * t45 + t6
     &3 * (t78 * t55 + t28)) + t56 * (t45 * (t47 * t79 + t64) + t67 * t5
     &7) * t12) - t9 * t1 * t70 * (t69 * t3 + ddreal(6) * t73 * t57)) + 
     &t2)
      t6 = ddreal(1) / ddreal(9)
      result = t6 * (t3 * t5 * ((t3 * (t51 * (t88 + t63 * (t87 * t55 + t
     &28)) + t56 * (t63 * (t48 * t84 * t55 + t46) + t88 * t47) * t12) + 
     &t9 * (t36 * (t85 * t49 + t71 * t84) + (t52 * t49 - t1) * (t50 * t8
     &4 + t85)) + t59 * (-t45 * (t82 - epinv * (t62 - t75 - t3) + ddreal
     &(12) * t20 + t81 * (ddreal(1) + t20) + ddreal(4)) - t87 * t66) + t
     &8 * (t55 * t73 * t84 * t29 + t45 * (-t72 * t83 + t80 + ddreal(6) *
     & t86 + ddreal(4))) + t74 * (t63 * (t16 + t9) - t75 - t3)) * struc7
     &MP + t4) + t2) / (ecossin)**2

           ampNonresonantHeavyImC7MP = result
       end function ampNonresonantHeavyImC7MP

       function ampNonresonantHeavyImC7PM()
           implicit none
           type(dd_complex) :: ampNonresonantHeavyImC7PM
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           type(dd_complex) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           type(dd_complex) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           type(dd_complex) ::  t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52
           type(dd_complex) ::  t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63
           type(dd_complex) ::  t64,t65,t66,t67,t68,t69,t7,t70,t71,t72,t73,t74
           type(dd_complex) ::  t75,t76,t77,t78,t79,t8,t80,t81,t82,t83,t84,t85
           type(dd_complex) ::  t86,t87,t88,t89,t9,t90,t91,t92,t93,t94,t95

           type(dd_complex) :: result

      t1 = intHLs250000x0122D4eps0()
      t2 = s25 + s26 + s56
      t3 = s16 + s26 - s34 + s56
      t4 = s13 - s45 - s46
      t5 = (gw)**2
      t6 = (t5)**2
      t7 = (gb)**2
      t8 = (t7)**2
      t9 = t3 * s25
      t10 = s12 + s16 + s26 - s34
      t11 = t3 * t2
      t12 = t11 * propW16
      t13 = t12 * t5 * t7
      t14 = t9 * propW16
      t15 = propZ25 * s25
      t16 = t15 * t6
      t17 = (s12 + s13 + s16 + s26 - s34 - s45 - s46) * t2
      t18 = ddreal(6)
      t19 = ddreal(3)
      t20 = s12 + s16 + s25
      t21 = -s16 + s25 + s34
      t22 = s12 + s16 + s26 - s34
      t23 = -s45 - s46
      t24 = (s12 + s13 + s16) * s25
      t25 = (-s13 + s25 + s26) * s34
      t26 = s13 * s16
      t27 = (t20 + s26) * s26 + s56 * t22 + t21 * t23 + t24 - t25 - t26
      t28 = ddreal(2) * s13
      t29 = s25 + s16 - s34
      t30 = ddreal(2) * s26
      t31 = s25 - s13 + s45 + s46
      t32 = t19 * s26
      t33 = s12 + s26
      t34 = ddreal(2) * t33
      t35 = t23 * (t30 + t29)
      t36 = t12 * (t34 + s25)
      t37 = (s56)**2
      t38 = -(s12 + s13 - s16) * s25 - (t28 - s25 + s12 - s16 - s26) * s
     &26 + (ddreal(2) * t31 - s12 + s16 - s34 + t32) * s56 - t25 - t26 -
     & t35 + t36 + ddreal(2) * t37
      t39 = t15 - ddreal(2)
      t40 = s12 + s13
      t41 = ddreal(2) * t40
      t42 = -t41 + s16
      t43 = s16 + s25
      t44 = ddreal(2) * s12
      t45 = t19 * s13
      t46 = -t44 - t45 + t43
      t47 = s16 - s34
      t48 = ddreal(2) * s25
      t49 = ddreal(4) * s26
      t31 = -t19 * t31 + t44 - t47 - t49
      t23 = t23 * (t32 + t48 + t47)
      t32 = t18 * t37 - ddreal(2) * (-t46 - s26) * s26 - ddreal(2) * t25
     & - ddreal(2) * t26 - ddreal(2) * t23 + ddreal(2) * t42 * s25 - ddr
     &eal(2) * t31 * s56 + t36 * t19
      t50 = t15 * t7
      t51 = t50 * t27
      t52 = t5 * t39
      t53 = t16 * t19 * t38
      t54 = intHLs250000x0112D2eps1()
      t55 = t19 * epinv
      t56 = t55 - ddreal(4)
      t25 = -ddreal(8) * t42 * s25 + ddreal(8) * (-t46 - s26) * s26 + dd
     &real(8) * t31 * s56 + ddreal(8) * t23 + ddreal(8) * t25 + ddreal(8
     &) * t26
      t31 = ddreal(12) * t15 * t38 - ddreal(12) * t36
      t36 = ddreal(24) * t37
      t46 = t55 * t39
      t57 = ddreal(9) * t16 * epinv
      t58 = intHLs250000x0111D2eps0()
      t59 = t19 * s25
      t60 = ddreal(4) * t47
      t61 = t18 * t33
      t62 = s13 - s45 - s46
      t63 = t19 * t62
      t64 = s25 + s26
      t65 = t47 * s13
      t66 = t64 * s34
      t67 = ddreal(4) * s12
      t68 = (s13 + s25 + s26) * s34
      t41 = ddreal(4) * (t41 + s16) * s25 + ddreal(4) * (t44 + t45 + t43
     & + s26) * s26 + ddreal(4) * (t63 + t44 - t21 - s56) * s56 + ddreal
     &(4) * t23 + ddreal(4) * t26 - ddreal(4) * t68 - ddreal(2) * t12 * 
     &(t61 + t60 + t59) + t15 * (-t18 * ((s13 + s16 + s25 + s26) * s26 -
     & t66) - ddreal(2) * (s16 * t19 + t41) * s25 - ddreal(2) * (t19 * (
     &s13 + s16 - s34 - s45 - s46) + s25 + t44 + t49 + s56) * s56 - ddre
     &al(2) * t23 - ddreal(2) * t65 - t67 * s26 + t12 * (ddreal(8) * t47
     & + t59 + t61))
      t49 = intHLs250000x0113D4eps0()
      t61 = ddreal(4) * t40
      t69 = t18 * t62
      t70 = ddreal(10)
      t71 = t15 * t19
      t33 = (ddreal(8) * s12 + ddreal(12) * s13) * s26 + t70 * ((-t43 - 
     &s26) * s26 + t66) - ddreal(18) * t37 + ddreal(2) * (-ddreal(5) * s
     &16 + t61) * s25 + ddreal(2) * (-ddreal(9) * s25 - ddreal(14) * s26
     & - ddreal(5) * t47 + t67 + t69) * s56 + ddreal(4) * t23 + ddreal(4
     &) * t65 + t71 * (t19 * t37 + t42 * s25 + (-ddreal(4) * s13 + t43 -
     & t44 + s26) * s26 + (t28 - t64) * s34 + (t47 + ddreal(4) * s26 - d
     &dreal(4) * s13 + ddreal(4) * s45 + ddreal(4) * s46 - t44 + t59) * 
     &s56 + t12 * (ddreal(4) * t33 + t47 + t48) - ddreal(2) * t26 - ddre
     &al(2) * t35) - t12 * t18 * (t34 + t21)
      t11 = ddreal(9) * t11 * t6 * (-propW16 * t47 + ddreal(1))
      t27 = ddreal(2) * t8 * t27
      t42 = t5 * t7
      t44 = intHLs250000x0122D4eps1()
      t66 = -ddreal(4) + epinv
      t62 = ddreal(2) * t62
      t70 = ddreal(2) * t22
      t72 = t70 + s25
      t73 = propW16 * t72 * t2
      t74 = intHLs250000x0111D2eps1()
      t75 = -ddreal(2) + epinv
      t34 = t12 * (t34 + t29)
      t76 = intHLs250000x0113D4eps1()
      t4 = (-t1 * (t18 * t9 * t5 * t7 * (propZ25 * t4 + (-ddreal(1) + (-
     &s12 - s16 - s26 + s34) * propZ25) * t2 * propW16) - ddreal(12) * t
     &13 * t10 - ddreal(9) * t16 * ((-t2 - t3) * s13 + s45 * t3 + (t2 + 
     &t3) * s46 + t2 * (t14 - s12 - s16 - s26 + s34 + s45)) - ddreal(18)
     & * t14 * t2 * propZ25 * t6 * t10 + ddreal(8) * t17 * t5 * t7 + ddr
     &eal(4) * t3 * t5 * t7 * t4 - t15 * t8 * (t3 * (-s13 + s45 + s46) +
     & t17) - t13 * t19 * propZ25 * (s25)**2) - t54 * (-t57 * t38 - t7 *
     & (t5 * (-t46 * t32 + t25 + t31 - t36) - t51 * t56)) + t44 * (t57 *
     & ((-t28 - t20 - s26) * s26 - (t62 + t22) * s56 - t24 - t26 - t35 +
     & t68 + t12 * t72) + t7 * (t5 * (epinv * (t3 * (t73 * t18 + t71 * (
     &t73 - t62)) - ddreal(4) * (ddreal(2) * t20 + t45) * s26 + ddreal(4
     &) * (s13 + ddreal(2) * t64) * s34 - ddreal(4) * (t63 + t70) * s56 
     &- ddreal(4) * t23 - ddreal(4) * t26 - ddreal(8) * (s26)**2 - ddrea
     &l(8) * t24) - t25 - t31 + t36) + t51 * t66))) * s25
      t9 = s25 * intHLs250000x0112D2eps0()
      t10 = s26 + s56
      t12 = t29 + ddreal(2) * t10
      t13 = t10 * t19 + t47 + t48
      t14 = t8 * t21
      t17 = t6 * t12
      t20 = ddreal(2) * t52 * t7
      t22 = -t20 * t13 + t15 * (t17 * t19 + t14)
      t24 = ddreal(2) * t47
      t29 = t18 * t10
      t45 = (propZ25 * t3 * t19 + ddreal(4)) * s25 + t24 + t29
      t17 = ddreal(9) * t17
      t48 = ddreal(2) * t42
      t12 = propZ25 * t12
      t10 = ddreal(12) * t10
      t57 = (-t12 * t18 + ddreal(8)) * s25
      t63 = t17 * epinv
      t68 = -t74 * t75 + t9
      t70 = -t76 * t75 - t49
      t72 = s25 * t7
      t73 = t19 * propW16
      t77 = t73 * t3
      t78 = -t77 + ddreal(2)
      t79 = ddreal(2) * propW16
      t80 = -t79 * t3 + ddreal(1)
      t81 = -t48 * t39
      t82 = -t80
      t83 = -t78
      t84 = t6 * t82
      t85 = t18 * propW16 * (ddreal(1) - t15)
      t86 = t15 + ddreal(2)
      t77 = t77 * t86 - ddreal(4)
      t87 = ddreal(9) * t84
      t88 = ddreal(12) * propW16
      t89 = t88 * t3
      t90 = t15 * t18
      t82 = t90 * t82
      t91 = t8 * t56
      t92 = t87 * epinv
      t93 = t8 * t66
      t79 = t79 * t2 + ddreal(1)
      t73 = t73 * t2 + ddreal(1)
      t94 = t6 * t79
      t95 = t15 * (-t94 * t19 + t8)
      t86 = t19 * (propW16 * t86 * t2 + t15) + ddreal(2)
      t94 = ddreal(9) * t94
      t79 = t90 * t79
      t88 = t88 * t2
      t90 = t94 * epinv
      t73 = t3 * (-t19 * (-t58 * (-t81 * t73 + t95) + t68 * (t20 * t73 +
     & t95)) + (t1 * (t48 * t86 + t15 * (t94 + t8)) + t44 * (t15 * (t90 
     &+ t93) + t48 * (epinv * t86 - t79 + t88 + ddreal(4))) - t54 * (t15
     & * (-t90 + t91) + t48 * (t46 * t73 - t79 + t88 + ddreal(4)))) * s2
     &5 - ddreal(8) * t72 * (t5 * (-t85 * t2 + t71 - ddreal(2)) + t50) *
     & t70) * struc7PM
      t1 = -(t19 * (-t58 * t22 + t22 * t68) - (t1 * (-t48 * t45 + t15 * 
     &(-t17 + t14)) + t44 * (t15 * (t14 * t66 - t63) - t48 * (epinv * t4
     &5 + t10 + t57 + t60)) - t54 * (t15 * (t14 * t56 + t63) - t48 * (t4
     &6 * t13 + t10 + t57 + t60))) * s25 + ddreal(8) * t72 * (t5 * ((-t1
     &2 * t19 + ddreal(4)) * s25 + t24 + t29) + t50 * t21) * t70) * stru
     &c2PM + (-t19 * (-t58 * (t81 * t78 + t15 * (t6 * t80 * t19 + t8)) +
     & t68 * (t20 * t83 + t15 * (-t84 * t19 + t8))) + (t1 * (t48 * t77 +
     & t15 * (t87 + t8)) + t44 * (t15 * (t92 + t93) + t48 * (epinv * t77
     & - t82 + t89 - ddreal(8))) - t54 * (t15 * (-t92 + t91) + t48 * (t4
     &6 * t83 - t82 + t89 - ddreal(8)))) * s25 - ddreal(8) * t72 * (t5 *
     & (-t85 * t3 - t71 + ddreal(4)) + t50) * t70) * t2 * (struc3PM + st
     &ruc5PM + struc8PM) + t73
      t3 = ddreal(1) / t3
      t2 = ddreal(1) / t2
      result = im * propW34 * (-(t19 * (-t58 * (-t53 + t7 * (t41 * t5 + 
     &t51)) - t74 * (t16 * (t18 * (-s25 * t40 - (t28 + s12) * s26 + (-t6
     &2 + s25 - s12 + s26 + s56) * s56 - t35 - t65 + t34) - t55 * t38) +
     & t7 * (t5 * (epinv * t41 - ddreal(2) * t39 * (t19 * (t34 + t37) - 
     &(t61 + s16) * s25 - (s13 * t18 + s26 + t43 + t67) * s26 + (t28 + t
     &64) * s34 + (-t69 - t67 + t59 + t30 - t47) * s56 - ddreal(2) * t26
     & - ddreal(2) * t23)) + t51 * t75)) + t9 * (-t53 + t7 * (t52 * t32 
     &+ t51))) - ddreal(4) * s25 * (t49 * (t42 * t33 + t15 * (t27 - t11)
     &) + t76 * (t15 * (-t11 * epinv + t27 * t75) + t42 * (epinv * t33 -
     & t25 - t31 + t36))) - t4) * struc1PM + ddreal(2) * t1) / (ecossin)
     &**2 / (gw)**2 / s25 * t2 * t3 / ddreal(18)

           ampNonresonantHeavyImC7PM = result
       end function ampNonresonantHeavyImC7PM

       function ampNonresonantHeavyImC7PP()
           implicit none
           type(dd_complex) :: ampNonresonantHeavyImC7PP

           type(dd_complex) :: result

      result = ddreal(0)

           ampNonresonantHeavyImC7PP = result
       end function ampNonresonantHeavyImC7PP

       function ampNonresonantHeavyMM()
           implicit none
           type(dd_complex) :: ampNonresonantHeavyMM
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           type(dd_complex) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           type(dd_complex) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           type(dd_complex) ::  t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52
           type(dd_complex) ::  t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63
           type(dd_complex) ::  t64,t65,t66,t7,t8,t9

           type(dd_complex) :: result

      t1 = intHLs250000x0112D2eps0()
      t2 = s25 + s26 + s56
      t3 = s16 + s26 - s34 + s56
      t4 = s26 + s56
      t5 = ddreal(2)
      t6 = (s25 + s26 + s56 - s13 + s45 + s46) * (s34 + s45 + s46 + s56)
      t7 = (s12 * t2 - t5 * s45 * t4 + s13 * t3 - s16 * (s45 + s46 + s56
     &) - (s25 - s34) * s45 - s46 * (s26 - s34) - (s26 - s34 + s46) * s5
     &6 - (s56)**2) * s12 + t6 * t3
      t8 = t5 * s16
      t9 = s16 + s25
      t10 = s12 - s13
      t11 = -s13 + s25 + s26
      t12 = s16 + s26
      t13 = s25 + s26
      t14 = t11 * t12
      t15 = -s16 + s25 + s34
      t16 = s25 + s26 - s13 + s34
      t17 = s16 + s26 - s34
      t18 = t16 * t17
      t19 = s12 * t15 + t18
      t20 = -s12 - s13 + s25 + s26 + s34
      t21 = t5 * s45
      t22 = t21 + t20
      t23 = ddreal(3)
      t24 = t23 * s12
      t25 = t5 * s12
      t26 = s25 - s13 + s16
      t27 = ddreal(5)
      t28 = t23 * s26
      t29 = s45 + s16
      t30 = (t29 * t5 - s12 - s13 + s25 - s34 + s46 + t28) * s46
      t31 = t5 * (s26 + s45 + s46)
      t32 = t24 - t31 - t26
      t6 = (-t21 - s25 - s56 + s13 - s46 + s16 - s34 + s12) * s12 + t6
      t21 = s13 * s16
      t33 = (s45)**2
      t34 = (s46)**2
      t35 = t11 * s34
      t36 = t17 * t33
      t37 = t17 * t34
      t38 = (-s25 * (t8 + s12) - s26 * (t5 * t9 + t10) + t21) * s12
      t10 = ((-t32 + s56) * s56 - s12 * t10 - (t24 + s13) * s16 - (t25 -
     & s16) * s25 + (-s12 * t27 + s26 + t26) * s26 + (t24 + t12 - s34) *
     & s34 + (t28 + t8 + s25 - s13 - s34 + s45) * s45 + t30) * s56
      t39 = t19 * s45
      t17 = t17 * t22 * s46
      t40 = propW16 * t2
      t41 = s12 * (s26)**2 - t40 * t3 * t6
      t42 = t5 * t41
      t43 = t42 - (-t35 - s12 * (-t13 * t5 + s13) + t14) * s34 - t36 - t
     &37 - t38 - t10 - t39 - t17
      t44 = s12 * (-t13 * t23 + s13) - t14
      t45 = t5 * s25
      t46 = -t45 - s26 + s16 - s34
      t47 = ddreal(4)
      t48 = t47 * s12
      t49 = ddreal(7)
      t50 = s12 * s13
      t51 = (s12)**2
      t8 = t5 * t51 + (t48 + s13) * s16 + (t24 - s16) * s25 - (-s12 * t4
     &9 + s26 + t26) * s26 - (t48 + t12 - s34) * s34 - (t28 + t8 + s25 +
     & s12 - s13 - s34 + s45) * s45 - t30 - t50
      t24 = t48 - t31 - t26
      t6 = -s12 * t11 + (t20 + s45) * s45 + (t22 + s46) * s46 + (-t5 * (
     &s12 - s45 - s46) + t16 + s56) * s56 + t35 + t40 * t6
      t9 = (-s25 * (s16 * t23 + t25) - s26 * (t23 * t9 - s13 + t25) + t2
     &1) * s12
      t16 = ((-t24 + s56) * s56 - t8) * s56
      t20 = (s12 * t46 - t18) * s45
      t21 = ddreal(12) * t41
      t22 = propZ25 * s25
      t26 = ddreal(6) * t22 * t3 * t6
      t28 = (gb)**2
      t30 = (t28)**2
      t31 = (gw)**2
      t48 = (t31)**2
      t52 = t22 * t28
      t53 = t52 * t7
      t54 = t31 * (-t47 * ((-t35 - t44) * s34 + t17 + t36 + t37 + t9 + t
     &16 - t20) + t21 - t26)
      t55 = ddreal(9)
      t56 = t22 * t55 * t48
      t57 = intHLs250000x0121D2eps0()
      t58 = ddreal(1) / s25
      t59 = t55 * t48
      t60 = intHLs250000x0111D2eps1()
      t61 = -t5 + epinv
      t14 = -t14 + t50
      t25 = (t25 * t13 - t14 - t35) * s34
      t10 = -t42 + t39 + t17 + t25 + t38 + t10 + t36 + t37
      t39 = s16 - s34
      t50 = s16 + s26 + s45 + s46
      t62 = (s34)**2
      t12 = t12 * s26
      t63 = (s16 + s26 - s45 - s46) * s34
      t64 = t50 * s25
      t50 = s13 * t50
      t29 = t29 * s46
      t65 = s16 * s45
      t66 = t23 * s26 * (s45 + s46)
      t19 = -t42 + t19 * s45 + ((-t32 + s56) * s56 + t5 * (t65 + t29) + 
     &(-s26 * t27 - t23 * t39 - s12 + s13 - s46 - t45) * s12 + t33 + t34
     & - t62 + t12 + t63 + t64 - t50 + t66) * s56 + t17 + t25 + t36 + t3
     &7 + t38
      t25 = t47 * t39
      t24 = (-t24 + s56) * s56
      t18 = (-s12 * t46 + t18) * s45
      t11 = (ddreal(6) * s12 * t13 - t14 * t5) * s34 + t5 * ((t24 - t5 *
     & (-t65 - t29 + t51) + t66 + s12 * (-s25 * t23 - s26 * t49 + s13 + 
     &s45 - s46 - t25) + t12 + t33 + t34 - t50 - t62 + t63 + t64) * s56 
     &- t11 * t62 + t17 + t36 + t37 + t9 + t18) - ddreal(6) * t41 + t22 
     &* t23 * t3 * t6
      t12 = epinv * t10
      t13 = t30 * t7
      t14 = intHLs250000x0112D2eps1()
      t27 = ddreal(1) - epinv
      t29 = t22 * t48
      t32 = t5 * t31
      t33 = t27 * t28 * (t32 * t11 + t53)
      t34 = intHLs250000x0122D4eps1()
      t38 = intHLs250000x0113D4eps1()
      t39 = intHLs250000x0122D4eps0()
      t35 = (-t35 - t44) * s34
      t8 = t56 * t10 + t28 * (t31 * (t47 * ((t24 - t8) * s56 + t17 + t18
     & + t36 + t37 + t9 + t35) - t21 + t26) + t53)
      t18 = intHLs250000x0111D2eps0()
      t21 = intHLs250000x0113D4eps0()
      t24 = ddreal(8)
      t26 = intHLs250000x0121D2eps1() * epinv
      t42 = s12 + s16 + s26 - s34 - s45
      t2 = ddreal(1) / t2
      t42 = ddreal(1) / t42
      t44 = ddreal(1) / t3
      t45 = ddreal(1) / s12
      t2 = ddreal(1) / (gw)**2 * t42 * t2 / (ecossin)**2
      t6 = t2 * mytree * propW34 * (-t5 * t38 * (t55 * t29 * (t12 - t19)
     & - t33) + epinv2 * (-t59 * propZ25 * t43 + t28 * (t28 * propZ25 * 
     &t7 - t54 * t58)) - (t26 + t57 + t1) * (t56 * t43 - t28 * (-t54 + t
     &53)) + t14 * (-t55 * t29 * (epinv * t43 + t10) - t33) - t34 * (t55
     & * t29 * (t12 + t43) - t33) - t39 * t8 - t21 * (t22 * (ddreal(18) 
     &* t48 * t10 + t13 * t5) + t28 * t31 * (-ddreal(24) * t41 + t24 * (
     &t17 - t20 + t35 + t16 + t36 + t37 + t9) + ddreal(12) * t22 * t3 * 
     &t6)) - (t60 * (t5 * t61 * t28 * t31 * t11 + t22 * (t48 * (t12 * t5
     &5 - ddreal(18) * t19) + t13 * t61)) + t8 * t18) * t58) * t45 * t44
     & / ddreal(9)
      result = t2 * (s12 - s34 - s45 - s46) * (-t5 * (-t27 * t38 + t21) 
     &- t14 * t27 + t27 * t34 + t58 * (-t60 * t61 + epinv2 - t18) + t1 +
     & t26 - t39 + t57) * propW34 * ((struc9Step5 * (s12 - s34 - s45 - s
     &46 - s56) * t45 + struc8Step5) * (t22 * (t59 * (t40 * t5 + ddreal(
     &1)) + t30) + t32 * t28 * (t23 * (t40 * (t22 + t5) + t22) + t5)) + 
     &(t56 * (t4 * t5 + s16 + s25 - s34) + t28 * (t31 * ((ddreal(6) * pr
     &opZ25 * t3 + t24) * s25 + t25 + ddreal(12) * t4) - t52 * t15)) * s
     &truc4Step5 * t44) / ddreal(18) + t6

           ampNonresonantHeavyMM = result
       end function ampNonresonantHeavyMM

       function ampNonresonantHeavyMP()
           implicit none
           type(dd_complex) :: ampNonresonantHeavyMP

           type(dd_complex) :: result

      result = ddreal(0)

           ampNonresonantHeavyMP = result
       end function ampNonresonantHeavyMP

       function ampNonresonantHeavyPM()
           implicit none
           type(dd_complex) :: ampNonresonantHeavyPM

           type(dd_complex) :: result

      result = ddreal(0)

           ampNonresonantHeavyPM = result
       end function ampNonresonantHeavyPM

       function ampNonresonantHeavyPP()
           implicit none
           type(dd_complex) :: ampNonresonantHeavyPP
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           type(dd_complex) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           type(dd_complex) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           type(dd_complex) ::  t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52
           type(dd_complex) ::  t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63
           type(dd_complex) ::  t64,t65,t66,t67,t68,t69,t7,t70,t71,t72,t73,t74
           type(dd_complex) ::  t75,t76,t77,t78,t79,t8,t80,t81,t82,t83,t84,t85
           type(dd_complex) ::  t86,t87,t88,t89,t9,t90,t91,t92,t93,t94,t95,t96
           type(dd_complex) ::  t97,t98

           type(dd_complex) :: result

      t1 = intHLs250000x0112D2eps0()
      t2 = ddreal(3)
      t3 = ddreal(2)
      t4 = t2 * s25
      t5 = t4 * propZ25
      t6 = (gb)**2
      t7 = (gw)**2
      t8 = t6 * propZ25
      t9 = t8 * s25
      t10 = t7 * (-t3 + t5)
      t11 = t10 + t9
      t12 = intHLs250000x0113D4eps0()
      t13 = intHLs250000x0121D2eps0()
      t14 = intHLs250000x0122D4eps0()
      t15 = intHLs250000x0112D2eps1()
      t16 = ddreal(1) - epinv
      t17 = intHLs250000x0113D4eps1()
      t18 = intHLs250000x0122D4eps1()
      t19 = intHLs250000x0111D2eps0()
      t20 = -s16 + s25 + s34
      t21 = s26 + s56
      t22 = t21 * t3 + s16 + s25 - s34
      t23 = s16 - s34
      t24 = ddreal(4)
      t25 = s16 + s26 - s34 + s56
      t26 = intHLs250000x0111D2eps1()
      t27 = -t3 + epinv
      t28 = ddreal(1) / s25
      t29 = ddreal(1) / t25
      t30 = t26 * t27
      t31 = epinv2 * t28
      t32 = (t15 - t18) * t16
      t33 = intHLs250000x0121D2eps1() * epinv
      t34 = t28 * (t30 + t19) * t29
      t35 = t17 * t16
      t36 = t35 - t12
      t37 = t3 * t11 * t36
      t38 = s25 + s26 + s56
      t39 = s26 * t3 + s12 - s13 + s16 + s25 - s34 + s46 + s56
      t40 = (s12 + s16 + s25 + s26 - s34 + s56) * t38
      t41 = t40 * t7
      t42 = propZ25 * s25
      t43 = -ddreal(1) + t42
      t44 = t25 * t39
      t45 = t44 * t7
      t46 = t45 * propW16 * t38 * t43
      t47 = t9 * (t40 + t44)
      t48 = t42 * t7
      t49 = t48 * (t40 - t44)
      t50 = ddreal(6) * t46
      t51 = (s25)**2
      t52 = epinv * propW16 * t38
      t53 = t38 * t51
      t54 = t45 * epinv
      t55 = ddreal(8)
      t56 = ddreal(12) * t46
      t57 = -t43
      t58 = ddreal(6) * propW16
      t59 = t58 * t57
      t60 = t7 * (-t59 * t38 - t3 + t5) + t9
      t61 = t6 * propW34
      t62 = t61 * t60
      t63 = ddreal(1) / t39
      t64 = t7 * propW16 * t38 * t25
      t43 = t64 * t43
      t65 = t25 + t38
      t66 = -t25 + t38
      t67 = t48 * t66
      t68 = t25 * t3
      t69 = -t2 * t67 + t9 * t65 + t7 * (t24 * t38 - t68) + ddreal(6) * 
     &t43
      t70 = -ddreal(1) + epinv
      t71 = t70 * t38
      t72 = ddreal(1) - epinv
      t73 = t7 * t38
      t74 = t7 * t25
      t75 = s26 * t38
      t76 = t9 * (t75 - t44)
      t77 = t48 * (t75 + t44)
      t68 = t7 * ((propZ25 * (t68 * propW16 * t38 - t20) * t2 + t24) * s
     &25 - t3 * (s16 - s34 - s26 - s56) - t58 * t38 * t25) + t9 * t22
      t78 = t18 * t16
      t79 = t78 * t68
      t80 = (t33 + t1 + t13 - t14) * t69
      t48 = t15 * (t2 * t48 * (t25 * t72 + t71) + t24 * t73 * t16 - t3 *
     & t74 * t16 - t9 * (t25 * t70 + t71) + ddreal(6) * t64 * (t42 * t72
     & + epinv - ddreal(1)))
      t64 = t35 * t68
      t68 = t12 * t69
      t69 = t74 * epinv2
      t71 = t73 * t24 * epinv2
      t59 = t59 * t25 - t24 + t5
      t72 = -t7 * t59 + t9
      t59 = -t7 * t59
      t74 = t59 + t9
      t81 = t61 * t74
      t82 = s25 - s34 + s56 - s13
      t83 = s16 + s12 + s46
      t84 = t2 * s26
      t58 = t7 * (t24 * t83 + t55 * t82 + ddreal(12) * s26 + t5 * (t3 * 
     &(t44 * propW16 + s13 - s25 + s34 - s56) - t83 - t84) - t58 * t44) 
     &+ t9 * (t3 * t82 + t83 + t84)
      t82 = t61 * t58
      t83 = t19 * t82
      t57 = t61 * propW16 * t57
      t85 = ddreal(1) / (gw)**2
      t86 = s12 - s13
      t87 = t3 * s13
      t88 = ddreal(5) * s25
      t89 = t3 * s16
      t90 = t3 * s25
      t91 = s13 - s46
      t92 = (t86 + s16) * s16
      t93 = (s56)**2
      t94 = (s26)**2
      t95 = ddreal(5) * s16
      t96 = t3 * s12
      t97 = t2 * s12
      t98 = ddreal(5) * s26
      t4 = t37 - (-t10 * t28 - t8) * epinv2 - t11 * (-t33 + t32 - t1 - t
     &13 + t14) + t34 * t63 * (t7 * ((s25 * t55 - propZ25 * (t24 * t93 +
     & t3 * t51 + ddreal(6) * t94 + (s16 * t2 + s12 - t87) * s25 + (-s13
     & * t2 + ddreal(7) * s25 + t95 + t96) * s26 + (-t86 - t4 - t89 - t9
     &8 + s34) * s34 + (t90 + t84 + t23) * s46 + (-t2 * t91 + t23 * t24 
     &+ ddreal(6) * s25 + ddreal(10) * s26 + t96) * s56 + t92) * t2) * s
     &25 + t3 * ((-s13 * t24 + t95 + t96) * s25 + (-ddreal(5) * s13 + dd
     &real(7) * s16 + ddreal(13) * s25 + t97) * s26 + (-ddreal(7) * s26 
     &- t86 - t88 - t89 + s34) * s34 + (s25 * t24 + t23 + t98) * s46 + (
     &ddreal(11) * s25 + ddreal(17) * s26 + ddreal(6) * t23 - ddreal(5) 
     &* t91 + t97) * s56 + t92) + ddreal(14) * t93 + ddreal(20) * t94) +
     & t9 * (t3 * (t93 + t94 + t51) + (s16 + s12 - t87) * s25 - (s16 + s
     &13 - t88) * s26 - (s25 - s26 - s12 + s13 - t89 + s34) * s34 + (-s1
     &6 + s34 + s26 + t90) * s46 + (t24 * (s25 + s26) - t91) * s56 - t92
     &))
      t10 = t15 * t16
      t51 = ddreal(1) / t38
      t5 = (struc30PP * (t2 * t67 * epinv2 * t39 - t3 * t39 * ((-t68 + t
     &64) * s25 - t69) - t19 * (-t2 * t77 + t7 * (t75 * t24 + t44 * t3) 
     &- t50 + t76) + t26 * (t2 * t77 * epinv + t24 * t7 * (-t75 * epinv 
     &+ t44) + t3 * (-t54 + t76) - t56 - t76 * epinv - ddreal(6) * t7 * 
     &(t75 * t42 + t44 * (t52 * (-t42 + ddreal(1)) + t42)) + t73 * t55 *
     & s26) - t39 * s25 * (t8 * epinv2 * t65 - t48 + t79 + t80) - ddreal
     &(6) * t46 * epinv2 - t71 * t39) * t51 * t63 + struc35PP * t74 * (t
     &3 * s25 * (t17 * t70 + t12) - s25 * (-t10 + t33 + t78 + t1 + t13 -
     & t14) - epinv2 + t19 + t30) - struc36PP * ((-t3 * (-t35 + t12) + t
     &1 - t10 + t13 - t14 + t33 + t78) * s25 + epinv2) * (t7 * (-t5 + t2
     &4) + t9)) * t29 * t28
      t10 = (t3 * t62 * t36 - t62 * (-t1 - t13 + t14) - t61 * (t60 * (-t
     &33 + t32 - t31) + t28 * (t19 * (-t2 * t49 + t24 * t41 - t45 * t3 +
     & t47 + t50) - t26 * (t2 * t49 * epinv - t24 * t7 * (t40 * epinv + 
     &t44) + t3 * (t8 * ((t38 * (s12 + s16 + s26 - s34 + s56) + t44) * s
     &25 + t53) + t54) + t41 * t55 - t47 * epinv + ddreal(6) * t7 * (((t
     &38 * (-s12 - s16 - s26 + s34 - s56) + t44 * (-t52 + ddreal(1))) * 
     &s25 - t53) * propZ25 + t52 * t44) + t56)) * t29 * t63)) * struc20P
     &P * t51
      t1 = t85 * (t29 * (struc31PP * (-t3 * t81 * t36 - (t1 + t13) * t61
     & * t72 + t14 * t81 + t83 * t28 * t63 - t61 * (epinv2 * (t59 * t28 
     &+ t8) + t33 * t72 - t32 * t74 - t30 * t58 * t63 * t28)) - struc34P
     &P * (-t2 * t9 * epinv2 * t7 * propW34 * t66 - t3 * t61 * ((t68 - t
     &64) * s25 + t69) - t19 * t25 * t62 + t42 * epinv2 * t65 * propW34 
     &* (t6)**2 + t61 * (s25 * (-t48 + t80 + t79) - t25 * t26 * t60 * t2
     &7) + ddreal(6) * t43 * t61 * epinv2 + t71 * t61) * t51 * t28) - t1
     &0 + t61 * (t51 * (struc10PP * (-t37 + t11 * (-t33 + t32 - t31 - t1
     & - t13 + t14) - t34 * (t7 * ((-propZ25 * t22 * t2 + t24) * s25 + t
     &23 * t3 + ddreal(6) * t21) + t9 * t20)) + struc33PP * t4) + t5)) +
     & (ddreal(6) * t57 * (t16 * (-t15 + t18) + t1 + t13 - t14 + t31 + t
     &33) + ddreal(12) * t57 * t36 + t85 * t28 * t29 * t63 * (t30 * t82 
     &+ t83)) * struc32PP
      result = t1 / (ecossin)**2 / ddreal(9)

           ampNonresonantHeavyPP = result
       end function ampNonresonantHeavyPP

       function ampNonresonantHeavyReC4MM()
           implicit none
           type(dd_complex) :: ampNonresonantHeavyReC4MM

           type(dd_complex) :: result

      result = ddreal(0)

           ampNonresonantHeavyReC4MM = result
       end function ampNonresonantHeavyReC4MM

       function ampNonresonantHeavyReC4MP()
           implicit none
           type(dd_complex) :: ampNonresonantHeavyReC4MP
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           type(dd_complex) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           type(dd_complex) ::  t31,t32,t33,t34,t35,t36,t37,t38,t4,t5,t6,t7
           type(dd_complex) ::  t8,t9

           type(dd_complex) :: result

      t1 = intHLs250000x0121D2eps0()
      t2 = intHLs250000x0121D2eps1() * epinv
      t3 = (t2 + t1) * s25 + epinv2
      t4 = (gb)**2
      t5 = propZ25 * s25
      t6 = t5 * (gw)**2
      t7 = s25 + s26 + s56
      t8 = s16 + s26 - s34 + s56
      t9 = ddreal(2) * propW16
      t10 = ddreal(6) * propW16
      t11 = ddreal(3) * t6
      t12 = intHLs250000x0112D2eps1()
      t13 = intHLs250000x0113D4eps1()
      t14 = intHLs250000x0122D4eps1()
      t15 = ddreal(2) * t13
      t16 = (-t15 + t2 + t1 + t12 - t14) * s25 + epinv2
      t17 = t9 * t8 - ddreal(1)
      t18 = s25 + s16 - s34
      t19 = s26 + s56
      t20 = t18 + ddreal(2) * t19
      t21 = (s34)**2
      t22 = (s16)**2
      t23 = (s25)**2 + t21 + t22
      t24 = (s56)**2
      t21 = s25 * (s16 * s25 - s25 * s34 + t21 + t22)
      t22 = s12 * t7
      t25 = s26 * propW16
      t26 = t25 * ((s25 + s56) * s16 + (-s25 - s56) * s34 + s56 * (s25 +
     & s26))
      t27 = (s26)**2
      t28 = s16 * s34
      t28 = (((-s16 + s34) * s56 + t28) * s25 - s26 * t27 + (t28 - t24) 
     &* s26 + t28 * s56) * propW16 - s45 * t19 - s46 * t19
      t29 = t27 * propW16 * (s25 + s16 - s34)
      t30 = t10 * t8 + t5 - ddreal(4)
      t31 = -ddreal(4) * s25 + ddreal(2) * s34 - ddreal(6) * t19 - s16 *
     & (t5 + ddreal(2)) + t5 * (s25 + s34)
      t32 = -ddreal(1) + t5
      t33 = t32 * s56
      t32 = t32 * s26
      t34 = s25 * (-ddreal(1) + (s25 + s45 + s46) * propZ25)
      t35 = -s45 - s46
      t36 = t4 * t30
      t37 = t11 * t17 + t36
      t38 = (t2 + t1 - t14) * t4
      t9 = (t4 * (t10 * t7 + t5 + ddreal(2)) + t11 * (t9 * t7 + ddreal(1
     &))) * t8
      t1 = t16 * (-struc2MP * (t11 * t20 - t4 * t31) + t7 * t37 * (struc
     &3MP + struc8MP)) - t7 * (epinv2 * (t4 * (-t5 + ddreal(4)) + t11) +
     & s25 * (ddreal(6) * t6 * t13 * t17 + ddreal(6) * (t4 * t8 + t6 * t
     &8) * t14 * propW16 - t12 * t37 - t38 * t5 + ddreal(3) * t6 * (t2 +
     & t1 - t14) + ddreal(4) * t38 + t36 * t15)) * struc5MP + t9 * (-t15
     & + t12 - t14) * s25 * struc7MP
      t2 = ddreal(1) / t7
      t10 = ddreal(1) / t8
      result = -ddreal(2) / ddreal(3) * propW34 * (t9 * t3 * (struc15MP 
     &+ struc41MP + struc57MP) - t16 * (t4 * (ddreal(12) * t26 + ddreal(
     &9) * t29 - ddreal(6) * t28 + ddreal(3) * (propW16 * t23 + s25) * s
     &56 + ddreal(3) * t24 * (propW16 * t18 + ddreal(1)) + ddreal(3) * t
     &25 * t23 + ddreal(3) * t21 * propW16 + s13 * t31 + (t34 + t33 + t3
     &2) * s16 + (s26 * (propZ25 * t19 - ddreal(1)) + t5 * (s26 - s45 - 
     &s46)) * s25 + (-t34 - t33 - t32) * s34 - t27 + t22 * t30 - ddreal(
     &2) * s16 * t35 + ddreal(2) * s26 * s56 + ddreal(2) * s34 * t35 - d
     &dreal(4) * s25 * t35) + t11 * (ddreal(4) * t26 + (s26 * t23 + s56 
     &* t23 + t18 * t24 + t21) * propW16 - s13 * t20 + s45 * t18 + s46 *
     & t18 + s56 * (s25 + s26) + t24 + t22 * t17 - ddreal(2) * t28 + ddr
     &eal(3) * t29)) * struc1MP - ddreal(2) * t1 + ddreal(6) * t8 * t7 *
     & (t6 + t4) * t3 * propW16 * (struc14MP + struc40MP + struc56MP)) /
     & (ecossin)**2 / s25 * t2 * t10

           ampNonresonantHeavyReC4MP = result
       end function ampNonresonantHeavyReC4MP

       function ampNonresonantHeavyReC4PM()
           implicit none
           type(dd_complex) :: ampNonresonantHeavyReC4PM
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           type(dd_complex) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           type(dd_complex) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           type(dd_complex) ::  t42,t43,t44,t5,t6,t7,t8,t9

           type(dd_complex) :: result

      t1 = s12 + s13
      t2 = ddreal(2)
      t3 = t2 * s13
      t4 = s25 + s16 - s34
      t5 = t2 * s26
      t6 = t5 + t4
      t7 = s25 + s26
      t8 = s13 - s45
      t9 = s25 + s26 + s56
      t10 = s16 + s26 - s34 + s56
      t11 = s13 * s16
      t12 = propW16 * (t2 * (s12 + s26) + t4) * t9 * t10
      t13 = t1 * s25
      t14 = t7 * s46
      t15 = intHLs250000x0121D2eps0()
      t16 = intHLs250000x0112D2eps1()
      t17 = s13 - s45 - s46
      t18 = -s16 + s34
      t19 = s45 + s46
      t20 = t18 * s13 - (s12 + t3) * s26 + (-t17 * t2 - s12 + s25 + s26 
     &+ s56) * s56 + t6 * t19 + t12 - t13
      t21 = intHLs250000x0113D4eps1()
      t22 = intHLs250000x0122D4eps1()
      t1 = s16 + ddreal(4) * t1
      t23 = ddreal(6)
      t24 = t23 * s13
      t25 = ddreal(4) * s12
      t26 = ddreal(3) * s26
      t27 = s25 * t2 - t18 + t26
      t28 = ddreal(3) * s25
      t29 = s25 - s16 + s34
      t30 = (s12 + s13 + s16) * s25
      t31 = (-s12 - s16 - s25 - s26) * s26 + (s25 - s13 + s26) * s34 - (
     &s12 + s16 + s26 - s34) * s56 + t19 * t29 + t11 - t30
      t32 = (t25 + t24 + s25 + s16 + s26) * s26
      t7 = (t3 + t7) * s34
      t5 = (t17 * t23 - t18 + t25 - t28 - t5) * s56
      t17 = (s56)**2
      t19 = t2 * (-t19 * t27 + t11)
      t33 = (propZ25 * t31 + t1) * s25 + t32 - t7 + t5 - ddreal(3) * t12
     & - ddreal(3) * t17 + t19
      t34 = propZ25 * s25
      t35 = intHLs250000x0121D2eps1() * epinv
      t36 = (t35 + t15) * s25 + epinv2
      t37 = (gb)**2
      t38 = (gw)**2
      t39 = t34 * t38
      t40 = t34 + t2
      t41 = ddreal(3) * t39
      t42 = s26 + s56
      t43 = t2 * propW16
      t44 = t23 * propW16
      t4 = ((-t2 * t21 + t15 + t16 - t22 + t35) * s25 + epinv2) * propW3
     &4 * ((t37 * (-s34 * t2 + t23 * t42 + ddreal(4) * s25 + s16 * t40 -
     & t34 * (s25 + s34)) + t41 * (t2 * t42 + t4)) * struc2PM - t10 * (t
     &37 * (t44 * t9 + t2 + t34) + t41 * (t43 * t9 + ddreal(1))) * struc
     &7PM + (t37 * (t44 * t10 + t34 - ddreal(4)) + t41 * (t43 * t10 - dd
     &real(1))) * t9 * (-struc3PM - struc5PM - struc8PM))
      t15 = (t37 * t40 + t41) * ((t35 + t15) * s25 + epinv2) * t10 * pro
     &pW34 * (struc27PM + struc43PM - struc57PM)
      t10 = ddreal(1) / t10
      t9 = ddreal(1) / t9
      t35 = ddreal(1) / ddreal(3)
      result = t35 * (-t2 * t15 - (s25 * ((ddreal(12) * t39 * t20 - ddre
     &al(4) * t37 * t33) * t21 * propW34 - t16 * t37 * (-t2 * propW34 * 
     &(t1 * s25 + t34 * t31 - ddreal(3) * t17 + t19 + t32 + t5 - t7) + t
     &12 * t23 * propW34)) + t2 * propW34 * t37 * (t36 * (t2 * (-s45 * t
     &27 + t11) + ddreal(3) * (s26)**2 - ddreal(3) * t12 - ddreal(3) * t
     &17 + (propZ25 * (-(s12 + s25) * s26 + (s25 - s13) * s34 + s45 * t2
     &9 - (s12 + s16 - s34 - s46) * s56 + t11 + t14 - t30) + t1) * s25 +
     & (ddreal(3) * s16 + s25 + t24 + t25) * s26 - (s25 + t26 + t3) * s3
     &4 + (t23 * t8 + ddreal(4) * s12 - ddreal(4) * s46 - t18 - t28) * s
     &56 - ddreal(4) * t14) - s25 * t22 * t33) + t23 * t34 * propW34 * t
     &38 * (t36 * ((t3 + s16 + s12 + s26) * s26 - (s13 + s26) * s34 - t6
     & * s45 + (t2 * t8 + s12 - s25 - s46 - s56) * s56 + t11 - t12 + t13
     & - t14) + s25 * (-t16 + t22) * t20)) * struc1PM - ddreal(4) * t4) 
     &/ (ecossin)**2 / s25 * t9 * t10

           ampNonresonantHeavyReC4PM = result
       end function ampNonresonantHeavyReC4PM

       function ampNonresonantHeavyReC4PP()
           implicit none
           type(dd_complex) :: ampNonresonantHeavyReC4PP

           type(dd_complex) :: result

      result = ddreal(0)

           ampNonresonantHeavyReC4PP = result
       end function ampNonresonantHeavyReC4PP

       function ampNonresonantHeavyReC7MM()
           implicit none
           type(dd_complex) :: ampNonresonantHeavyReC7MM

           type(dd_complex) :: result

      result = ddreal(0)

           ampNonresonantHeavyReC7MM = result
       end function ampNonresonantHeavyReC7MM

       function ampNonresonantHeavyReC7MP()
           implicit none
           type(dd_complex) :: ampNonresonantHeavyReC7MP
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           type(dd_complex) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           type(dd_complex) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           type(dd_complex) ::  t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52
           type(dd_complex) ::  t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63
           type(dd_complex) ::  t64,t65,t66,t67,t68,t69,t7,t70,t71,t72,t73,t74
           type(dd_complex) ::  t75,t76,t77,t78,t79,t8,t80,t81,t82,t83,t84,t85
           type(dd_complex) ::  t86,t87,t88,t89,t9,t90,t91,t92,t93,t94

           type(dd_complex) :: result

      t1 = intHLs250000x0112D2eps0()
      t2 = s25 + s16 - s12
      t3 = ddreal(2)
      t4 = t3 * s13
      t5 = s34 - s25 - s16
      t6 = t3 * s26
      t7 = s34 - s16 + s12
      t8 = s45 + s46 + s25 - s13
      t9 = ddreal(3)
      t10 = t9 * s26
      t11 = s12 + s26
      t12 = t3 * t11
      t13 = s25 + s26 + s56
      t14 = s34 - s56 - s16 - s26
      t15 = -s45 - s46
      t16 = s13 * s16
      t17 = t15 * (t6 - t5)
      t18 = (s12 + s13 - s16) * s25
      t19 = (-s13 + s25 + s26) * s34
      t20 = propW16 * t13
      t21 = t20 * (t12 + s25)
      t22 = t21 * t14
      t23 = (s56)**2
      t24 = t3 * t23 - (t4 - t2 - s26) * s26 + (t3 * t8 + t10 - t7) * s5
     &6 - t16 - t17 - t18 - t19 - t22
      t25 = s12 + s13
      t26 = t3 * t25
      t27 = t26 - s16
      t28 = s16 + s25
      t29 = t9 * s13
      t30 = t3 * s12
      t31 = t30 + t29 - t28
      t32 = s34 - s16
      t33 = t3 * s25
      t34 = ddreal(4) * s26
      t35 = s12 + s16 + s25
      t36 = -s16 + s25 + s34
      t37 = s34 - s16 - s12 - s26
      t38 = t15 * t36
      t39 = (s12 + s13 + s16) * s25
      t40 = (-t35 - s26) * s26 + s56 * t37 + t16 + t19 - t38 - t39
      t41 = (s26)**2
      t42 = t15 * (t33 + t10 - t32)
      t8 = (-t8 * t9 + t30 + t32 - t34) * s56
      t43 = propZ25 * t40
      t44 = ddreal(6)
      t45 = (gb)**2
      t46 = (t45)**2
      t47 = (gw)**2
      t48 = (t47)**2
      t49 = t9 * s25
      t50 = t49 * t47 * propZ25
      t51 = intHLs250000x0113D4eps0()
      t52 = ddreal(4) * t25
      t53 = s25 + s26
      t54 = s45 + s46 - s13
      t55 = ddreal(9)
      t56 = t44 * t54
      t57 = ddreal(4) * s12
      t58 = t3 * t54
      t59 = t20 * (-t3 * t37 + s25)
      t60 = t32 * s13
      t61 = ddreal(10)
      t62 = t49 * propZ25
      t61 = (-ddreal(8) * s12 - ddreal(12) * s13) * s26 + t14 * (-t20 * 
     &(t12 + t36) * t44 - t62 * (t59 + t58)) + ddreal(18) * t23 - t3 * (
     &(-ddreal(5) * s16 + t52) * s25 + (-s25 * t55 - ddreal(14) * s26 + 
     &ddreal(5) * t32 - t56 + t57) * s56) + t61 * ((t28 + s26) * s26 - s
     &34 * t53) - ddreal(4) * t42 + ddreal(4) * t60
      t63 = propZ25 * s25
      t64 = t63 * t45
      t65 = t64 * t40
      t66 = t63 * t55
      t67 = t66 * t48
      t68 = t67 * t24
      t69 = intHLs250000x0122D4eps0()
      t70 = t44 * s13
      t59 = t59 * t14
      t71 = (s13 + s25 + s26) * s34
      t72 = (t4 + t35 + s26) * s26 + (-t58 - t37) * s56 + t16 + t17 + t3
     &9 + t59 - t71
      t35 = t3 * (-(t3 * t53 + s13) * s34 + t16 + t42) + t9 * (-t63 * t7
     &2 + t59) + (ddreal(4) * t35 + t70) * s26 + (-ddreal(4) * t37 - t56
     &) * s56 + ddreal(4) * t39 + ddreal(4) * t41
      t37 = intHLs250000x0111D2eps0()
      t39 = ddreal(4) * t32
      t2 = -t3 * ((t26 + s16) * s25 + (t30 + t29 + t28 + s26) * s26 + (-
     &t54 * t9 - s56 + t30 - t36) * s56 + t16 + t42 - t71) - t20 * (t11 
     &* t44 - t39 + t49) * t14 + t63 * (t3 * (t20 * t32 * t14 - t23) - (
     &t2 + s26) * s26 + (-t33 - t10 + t7) * s56 - t16 + t18 + t38 + t71)
      t7 = intHLs250000x0112D2eps1()
      t10 = t3 * (s45 + s46 + s56 + s25 - s13 + s26) + t21
      t11 = s56 + s26
      t18 = s25 + s26 + s45 + s46 + s56
      t26 = s25 * t27 + (t31 - s26) * s26 + t16 + t19 + t42 + t8
      t29 = (t62 * t10 + t21 * t44) * t14
      t30 = -ddreal(1) + epinv
      t38 = t63 * t48
      t54 = intHLs250000x0113D4eps1()
      t59 = -t3 + epinv
      t71 = intHLs250000x0122D4eps1()
      t73 = -t30
      t74 = intHLs250000x0111D2eps1()
      t12 = t20 * (t12 - t5) * t14
      t75 = t9 * epinv
      t76 = t63 * t47
      t2 = (t37 * propW34 * (t2 * t45 - t50 * t72) - s25 * t1 * propW34 
     &* (-t45 * (-t23 * t44 + t3 * (t8 + t19 + t42 + t16 - t41) + t3 * t
     &31 * s26 + (t3 * t27 + t43) * s25 + t22 * t9) + t50 * t24) + propW
     &34 * t74 * (t45 * (epinv * t2 + t3 * ((t43 + t52 + s16) * s25 + (t
     &57 + t70 + t28 + s26) * s26 - (t4 + t53) * s34 + (-t6 - t49 + t57 
     &- t56 - t32) * s56) + t44 * (t12 - t23) + ddreal(4) * t16 + ddreal
     &(4) * t42) + t76 * (t44 * (s25 * t25 + (t4 + s12) * s26 + (-t58 - 
     &s25 + s12 - s26 - s56) * s56 + t17 - t60 + t12) - t75 * t72))) * t
     &47
      t4 = t11 * t44
      t6 = t3 * t32
      t8 = (propZ25 * t36 - ddreal(4)) * s25 - t4 + t6
      t5 = t3 * t11 - t5
      t12 = ddreal(1) / (gw)**2
      t16 = t45 * t12
      t17 = t5 * t9
      t19 = -t45 * t8 + t50 * t5
      t25 = propZ25 * t14
      t27 = t46 * t12
      t28 = t47 * t5
      t31 = t28 * t55
      t41 = t3 * ((-t25 * t9 + ddreal(4)) * s25 + t4 - t6)
      t42 = t64 * t36
      t25 = (t25 * t44 - ddreal(8)) * s25
      t43 = ddreal(12) * t11
      t49 = t42 * t12
      t52 = ddreal(1) / s25
      t56 = t74 * t59
      t57 = t45 * t37
      t58 = t54 * t59
      t60 = t69 * t45
      t70 = ddreal(1) - t63
      t72 = t44 * propW16
      t74 = t72 * t70 * t14
      t77 = -t3 * propW16 * t14 - ddreal(1)
      t72 = t72 * t14
      t78 = -t72 + t63 - ddreal(4)
      t79 = t45 * t78
      t80 = t63 + t3
      t81 = t47 * t77
      t82 = t3 * (-t9 * propW16 * t80 * t14 - ddreal(4))
      t83 = t64 * t12
      t84 = t63 * t44
      t85 = t83 * t73
      t86 = t20 * t3 + ddreal(1)
      t87 = t20 * t44
      t88 = t87 + t63 + t3
      t70 = t87 * t70
      t87 = t45 * t88 + t50 * t86
      t89 = t20 * t80 + t63
      t90 = t47 * t86
      t91 = t90 * t55
      t92 = t3 * t45 * (t89 * t9 + t3)
      t93 = ddreal(1) / t13
      t94 = ddreal(1) / t14
      t2 = (t3 * propW34 * s25 * (t51 * (t68 - t45 * (-t47 * t61 + t65))
     & + t54 * (t38 * t24 * (epinv * t55 - ddreal(18)) + t45 * (t47 * (e
     &pinv * t61 + t14 * (t63 * t10 * t44 + ddreal(12) * t21) - ddreal(2
     &4) * t23 + ddreal(8) * t26) - t65 * t59))) + t9 * t2 + (t69 * t45 
     &* propW34 * (t35 * t47 - t65) + t7 * (-t55 * t38 * propW34 * t24 *
     & t30 + t45 * (t47 * (epinv * (t44 * (-t23 * t9 + t26) + t62 * t40 
     &+ t22 * t55) * propW34 + propW34 * (ddreal(4) * t3 * (-s12 * t13 -
     & s25 * t15) + ddreal(4) * t9 * (-s26 * t15 + (s25 + s45 + s46) * s
     &56 + t23) - ddreal(4) * s13 * (t11 * t9 - t32 + t33) + ddreal(4) *
     & t18 * s16 + ddreal(4) * s26 * t53 - ddreal(4) * t18 * s34 + ddrea
     &l(4) * t34 * s56 - t29)) - t65 * propW34)) + propW34 * t71 * (-t68
     & + t45 * (t47 * (epinv * t35 - ddreal(12) * t23 + ddreal(4) * t26 
     &+ t29) + t65 * t73))) * s25) * struc1MP * t12 * t52 * t93
      t4 = (-t3 * (-t51 * (-t41 * t45 + t63 * (t27 * t36 - t31)) + t58 *
     & (t67 * t5 - t45 * (-t41 * t47 + t42)) * t12) + t9 * (t1 * t19 + t
     &52 * (-t56 * t19 + t57 * t8)) + t7 * (t55 * t76 * t5 * t30 + t45 *
     & (-t75 * t8 + t25 + t39 - t43 + t49)) + t71 * (t28 * t66 - t45 * (
     &t49 * t73 + epinv * ((t17 * propZ25 - ddreal(4)) * s25 - t4 + t6) 
     &+ t25 + t39 - t43)) + t60 * (s25 * ((t16 * t36 - t17) * propZ25 + 
     &ddreal(4)) + t4 - t6) - t31 * t37 * propZ25) * struc2MP * t94
      t4 = propW34 * t93 * ((t3 * (t51 * (t92 + t63 * (t91 + t27)) + t58
     & * (t92 * t47 + t63 * (t48 * t86 * t55 + t46)) * t12) + t9 * (-t1 
     &* t87 + t52 * (t56 * t87 + t57 * t88)) + t7 * (-t55 * t76 * t86 * 
     &t30 + t45 * (t44 * t89 - t75 * t88 + t83 + ddreal(4))) + t71 * (-t
     &90 * t66 + t45 * (-t85 + epinv * (-t70 + t62 - t3) - t84 * (ddreal
     &(1) + t20) - ddreal(12) * t20 - ddreal(4))) + t60 * (t63 * (t16 + 
     &t9) - t70 - t3) + t91 * t37 * propZ25) * struc7MP - t4)
      t5 = ddreal(1) / ddreal(9)
      result = t5 * (t3 * t4 - t94 * ((struc3MP + struc5MP + struc8MP) *
     & (t3 * (t7 * (-t55 * t76 * t77 * t30 + t45 * (-t72 * t80 - t75 * t
     &78 + t83 - ddreal(8))) + t71 * (-t81 * t66 + t45 * (-t85 + epinv *
     & (t74 - t62 + ddreal(4)) + propW16 * t14 * (t84 + ddreal(12)) + dd
     &real(8))) + t60 * (t63 * (t16 - t9) + t74 + ddreal(4))) * propW34 
     &+ t44 * propW34 * (t50 * t77 + t79) * t52 * (t56 + t37) + ddreal(4
     &) * propW34 * (t51 * (t82 * t45 + t63 * (t81 * t55 + t27)) + t58 *
     & (t67 * t77 + t45 * (t47 * t82 + t64)) * t12) - t9 * t1 * propW34 
     &* (t76 * t44 * t77 + t79 * t3)) + t2)) / (ecossin)**2

           ampNonresonantHeavyReC7MP = result
       end function ampNonresonantHeavyReC7MP

       function ampNonresonantHeavyReC7PM()
           implicit none
           type(dd_complex) :: ampNonresonantHeavyReC7PM
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           type(dd_complex) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           type(dd_complex) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           type(dd_complex) ::  t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52
           type(dd_complex) ::  t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63
           type(dd_complex) ::  t64,t65,t66,t67,t68,t69,t7,t70,t71,t72,t73,t74
           type(dd_complex) ::  t75,t76,t77,t78,t79,t8,t80,t81,t82,t83,t84,t85
           type(dd_complex) ::  t86,t87,t88,t89,t9,t90,t91,t92,t93,t94,t95,t96

           type(dd_complex) :: result

      t1 = intHLs250000x0122D4eps0()
      t2 = s25 + s26 + s56
      t3 = s16 + s26 - s34 + s56
      t4 = s13 - s45 - s46
      t5 = (gb)**2
      t6 = (t5)**2
      t7 = (gw)**2
      t8 = (t7)**2
      t9 = propW16 * t2
      t10 = t5 * t7
      t11 = t10 * s25 * t3
      t12 = s12 + s16 + s26 - s34
      t13 = t10 * t9 * t3
      t14 = t2 + t3
      t15 = -t2 - t3
      t16 = propW16 * s25 * t3
      t17 = propZ25 * s25
      t18 = t17 * t8
      t19 = -s13 + s45 + s46
      t20 = t19 * t3
      t21 = (s12 + s13 + s16 + s26 - s34 - s45 - s46) * t2
      t22 = t17 * t6 * (t21 + t20)
      t23 = t18 * t9 * t3 * t12
      t24 = t10 * t2
      t25 = ddreal(3)
      t26 = (s25)**2
      t27 = t13 * t25 * propZ25 * t26
      t28 = intHLs250000x0122D4eps1()
      t29 = propZ25 * t14
      t26 = t3 * t26
      t30 = t9 * t3
      t31 = t18 * epinv
      t32 = -s12 - s13 + s45 + s46 + s56
      t19 = -ddreal(24) * t13 * (t17 * (-s12 - s26) + s12 + s26) - ddrea
     &l(16) * t24 * t32 - ddreal(12) * t10 * ((t3 * (-propZ25 * s26 + t9
     &) - t29 * s45 + t29 * s13 - t29 * s56 - t29 * s46 + propZ25 * s12 
     &* t2) * s25 - t26 * propZ25 * (t9 + ddreal(1)) + t30 * t12 * epinv
     &) - t22 * epinv + ddreal(9) * t31 * (s13 * t3 + s45 * t15 + s46 * 
     &t15 + t2 * (-t16 + s12 + s13 + s16 + s26 - s34)) - ddreal(18) * t2
     &3 * epinv - ddreal(4) * t5 * (t20 * t7 * epinv + ((t2 * t32 + t3 *
     & (s13 - s26 - s45 - s46 - s56)) * s25 - t26) * propZ25 * t5) + ddr
     &eal(8) * t10 * (t3 * (s13 - s25 - s26 - s45 - s46 - s56) + t21 * e
     &pinv) - ddreal(6) * t11 * epinv * (propZ25 * t19 + t9 * (propZ25 *
     & t12 + ddreal(1))) - t27 * epinv
      t20 = -s16 + s25 + s34
      t21 = -s45 - s46
      t26 = (-s13 + s25 + s26) * s34
      t29 = s13 * s16
      t32 = (s12 + s13 + s16) * s25 + (s12 + s16 + s25 + s26) * s26 + (s
     &12 + s16 + s26 - s34) * s56 + t21 * t20 - t26 - t29
      t33 = ddreal(2) * s13
      t34 = s25 + s16 - s34
      t35 = ddreal(2) * s26
      t36 = s25 - s13 + s45 + s46
      t37 = t25 * s26
      t38 = s12 + s26
      t39 = ddreal(2) * t38
      t40 = t21 * (t34 + t35)
      t41 = t30 * (s25 + t39)
      t42 = (s56)**2
      t43 = -(s12 + s13 - s16) * s25 - (-s25 + s12 - s16 + t33 - s26) * 
     &s26 + (ddreal(2) * t36 - s12 + s16 - s34 + t37) * s56 - t26 - t29 
     &- t40 + t41 + ddreal(2) * t42
      t44 = -ddreal(2) + t17
      t45 = s12 + s13
      t46 = ddreal(2) * t45
      t47 = s16 - t46
      t48 = s16 + s25
      t49 = ddreal(2) * s12
      t50 = t25 * s13
      t51 = t48 - t49 - t50
      t52 = s16 - s34
      t53 = ddreal(2) * s25
      t54 = ddreal(4) * s26
      t36 = -t25 * t36 + t49 - t52 - t54
      t21 = t21 * (t52 + t53 + t37)
      t37 = -ddreal(2) * (-t51 - s26) * s26 - ddreal(2) * t26 - ddreal(2
     &) * t29 - ddreal(2) * t21 + ddreal(2) * t47 * s25 - ddreal(2) * t3
     &6 * s56 + ddreal(6) * t42 + t41 * t25
      t55 = t17 * t5
      t56 = t55 * t32
      t57 = t18 * t25 * t43
      t58 = intHLs250000x0112D2eps1()
      t59 = t25 * epinv
      t60 = -ddreal(4) + t59
      t26 = -ddreal(8) * t47 * s25 + ddreal(8) * (-t51 - s26) * s26 + dd
     &real(8) * t36 * s56 + ddreal(8) * t21 + ddreal(8) * t26 + ddreal(8
     &) * t29
      t36 = ddreal(12) * t17 * t43 - ddreal(12) * t41
      t41 = ddreal(24) * t42
      t51 = t59 * t44
      t61 = intHLs250000x0111D2eps0()
      t62 = t25 * s25
      t63 = ddreal(4) * t52
      t64 = ddreal(6) * t38
      t65 = s13 - s45 - s46
      t66 = s25 + s26
      t67 = t52 * s13
      t68 = t66 * s34
      t69 = ddreal(4) * s12
      t46 = -ddreal(4) * (s16 + t46) * s25 - ddreal(4) * (t48 + t49 + t5
     &0 + s26) * s26 + ddreal(4) * (s13 + s25 + s26) * s34 - ddreal(4) *
     & (t25 * t65 - s56 - t20 + t49) * s56 - ddreal(4) * t21 - ddreal(4)
     & * t29 - t17 * (-ddreal(2) * (s16 * t25 + t46) * s25 - ddreal(2) *
     & (t25 * (s13 + s16 - s34 - s45 - s46) + s25 + t49 + t54 + s56) * s
     &56 - ddreal(2) * t21 - ddreal(2) * t67 - ddreal(6) * (s13 + s16 + 
     &s25 + s26) * s26 + ddreal(6) * t68 + t30 * (ddreal(8) * t52 + t62 
     &+ t64) - t69 * s26) + ddreal(2) * t30 * (t62 + t63 + t64)
      t50 = intHLs250000x0113D4eps0()
      t54 = ddreal(4) * t45
      t64 = ddreal(6) * t65
      t70 = -ddreal(10)
      t71 = t17 * t25
      t38 = -(ddreal(8) * s12 + ddreal(12) * s13) * s26 + t70 * ((-t48 -
     & s26) * s26 + t68) + ddreal(18) * t42 - ddreal(2) * (-ddreal(5) * 
     &s16 + t54) * s25 - ddreal(2) * (-ddreal(9) * s25 - ddreal(14) * s2
     &6 - ddreal(5) * t52 + t69 + t64) * s56 - ddreal(4) * t67 - ddreal(
     &4) * t21 + ddreal(6) * t30 * (t20 + t39) - t71 * (t25 * t42 + t47 
     &* s25 + (-ddreal(4) * s13 + t48 - t49 + s26) * s26 + (-t66 + t33) 
     &* s34 + (t52 + ddreal(4) * s26 - ddreal(4) * s13 + ddreal(4) * s45
     & + ddreal(4) * s46 - t49 + t62) * s56 + t30 * (ddreal(4) * t38 + t
     &52 + t53) - ddreal(2) * t40 - ddreal(2) * t29)
      t47 = ddreal(9) * t8 * t2 * t3 * (-propW16 * t52 + ddreal(1))
      t32 = ddreal(2) * t6 * t32
      t49 = intHLs250000x0111D2eps1()
      t68 = -ddreal(2) + epinv
      t30 = t30 * (t34 + t39)
      t39 = ddreal(2) * t44
      t70 = intHLs250000x0113D4eps1()
      t72 = s25 * intHLs250000x0112D2eps0()
      t73 = s25 * propW34
      t4 = t73 * (t1 * (ddreal(6) * t11 * (propZ25 * t4 + t9 * (-ddreal(
     &1) + (-s12 - s16 - s26 + s34) * propZ25)) - ddreal(12) * t13 * t12
     & - ddreal(9) * t18 * (s13 * t15 + s45 * t3 + s46 * t14 + t2 * (t16
     & - s12 - s16 - s26 + s34 + s45)) - ddreal(18) * t23 - ddreal(8) * 
     &t24 * (-s12 - s13 - s16 - s26 + s34 + s45 + s46) + ddreal(4) * t10
     & * t3 * t4 - t22 - t27) + t28 * t19 + t58 * (t5 * (t7 * (t51 * t37
     & - t26 - t36 + t41) + t56 * t60) - ddreal(9) * t31 * t43))
      t11 = s26 + s56
      t12 = t34 + ddreal(2) * t11
      t13 = t11 * t25 + t52 + t53
      t14 = t6 * t20
      t15 = t8 * t12
      t16 = t39 * t10
      t19 = t17 * (t15 * t25 + t14) - t16 * t13
      t22 = ddreal(2) * t52
      t23 = ddreal(6) * t11
      t24 = (propZ25 * t3 * t25 + ddreal(4)) * s25 + t22 + t23
      t15 = ddreal(9) * t15
      t27 = ddreal(2) * t10
      t12 = propZ25 * t12
      t11 = ddreal(12) * t11
      t31 = (ddreal(8) - ddreal(6) * t12) * s25
      t34 = t15 * epinv
      t53 = -ddreal(4) + epinv
      t74 = -t49 * t68 + t72
      t75 = -t70 * t68 - t50
      t76 = t5 * s25
      t77 = t25 * propW16 * t3
      t78 = ddreal(2) - t77
      t79 = -ddreal(2) * propW16 * t3 + ddreal(1)
      t80 = -t27 * t44
      t81 = -t79
      t82 = -t78
      t83 = t8 * t81
      t84 = ddreal(1) - t17
      t85 = ddreal(2) + t17
      t77 = t77 * t85 - ddreal(4)
      t86 = ddreal(9) * t83
      t87 = ddreal(6) * t17
      t81 = t87 * t81
      t88 = ddreal(12) * propW16 * t3
      t89 = t86 * epinv
      t90 = t6 * t60
      t91 = t6 * t53
      t92 = ddreal(1) + ddreal(2) * t9
      t93 = t9 * t25 + ddreal(1)
      t94 = t8 * t92
      t95 = t17 * (-t94 * t25 + t6)
      t85 = t25 * (t9 * t85 + t17) + ddreal(2)
      t94 = ddreal(9) * t94
      t87 = t87 * t92
      t92 = ddreal(12) * t9
      t96 = t94 * epinv
      t9 = t3 * (t25 * (t61 * (-t80 * t93 + t95) - t74 * (t16 * t93 + t9
     &5)) - ddreal(8) * t76 * (t7 * (-ddreal(6) * t9 * t84 + t71 - ddrea
     &l(2)) + t55) * t75 + (t1 * (t17 * (t94 + t6) + t27 * t85) + t28 * 
     &(t27 * (epinv * t85 - t87 + t92 + ddreal(4)) + t17 * (t96 + t91)) 
     &- t58 * (t27 * (t51 * t93 - t87 + t92 + ddreal(4)) + t17 * (-t96 +
     & t90))) * s25) * struc7PM
      t1 = propW34 * ((t25 * (-t61 * t19 + t19 * t74) + ddreal(8) * t76 
     &* (t7 * ((-t12 * t25 + ddreal(4)) * s25 + t22 + t23) + t55 * t20) 
     &* t75 - (t1 * (t17 * (-t15 + t14) - t27 * t24) + t28 * (-t27 * (ep
     &inv * t24 + t11 + t31 + t63) + t17 * (t14 * t53 - t34)) - t58 * (t
     &17 * (t14 * t60 + t34) - t27 * (t51 * t13 + t11 + t31 + t63))) * s
     &25) * struc2PM + (t25 * (t61 * (t17 * (t8 * t79 * t25 + t6) + t80 
     &* t78) - t74 * (t17 * (-t83 * t25 + t6) + t16 * t82)) - ddreal(8) 
     &* t76 * (t7 * (-ddreal(6) * propW16 * t84 * t3 - t71 + ddreal(4)) 
     &+ t55) * t75 + (t1 * (t17 * (t6 + t86) + t27 * t77) + t28 * (t17 *
     & (t89 + t91) + t27 * (epinv * t77 - t81 + t88 - ddreal(8))) - t58 
     &* (t17 * (-t89 + t90) + t27 * (t51 * t82 - t81 + t88 - ddreal(8)))
     &) * s25) * t2 * (-struc3PM - struc5PM - struc8PM) - t9)
      t3 = ddreal(1) / t3
      t2 = ddreal(1) / t2
      result = ((-t25 * (-(t61 * (t5 * (t46 * t7 - t56) + t57) + t72 * (
     &t5 * (t7 * t44 * t37 + t56) - t57)) * propW34 + t49 * (t5 * (-t7 *
     & propW34 * (epinv * t46 + t39 * (t25 * (t42 + t30) - (s16 + t54) *
     & s25 - (ddreal(6) * s13 + t48 + t69 + s26) * s26 + (t66 + t33) * s
     &34 + (-t52 + t35 + t62 - t69 - t64) * s56 - ddreal(2) * t21 - ddre
     &al(2) * t29)) + t56 * t68 * propW34) - t18 * propW34 * (ddreal(6) 
     &* s25 * t45 + ddreal(6) * (s12 + t33) * s26 + ddreal(6) * (ddreal(
     &2) * t65 - s25 + s12 - s26 - s56) * s56 + ddreal(6) * t40 + ddreal
     &(6) * t67 - ddreal(6) * t30 + t59 * t43))) - ddreal(4) * t73 * (t5
     &0 * (-t10 * t38 + t17 * (-t47 + t32)) + t70 * (t17 * (-t47 * epinv
     & + t32 * t68) - t10 * (epinv * t38 + t26 + t36 - t41))) + t4) * st
     &ruc1PM / ddreal(18) + t1 / ddreal(9)) / (ecossin)**2 / (gw)**2 / s
     &25 * t2 * t3

           ampNonresonantHeavyReC7PM = result
       end function ampNonresonantHeavyReC7PM

       function ampNonresonantHeavyReC7PP()
           implicit none
           type(dd_complex) :: ampNonresonantHeavyReC7PP

           type(dd_complex) :: result

      result = ddreal(0)

           ampNonresonantHeavyReC7PP = result
       end function ampNonresonantHeavyReC7PP

       function ampNonresonantLightFullImC4MM()
           implicit none
           type(dd_complex) :: ampNonresonantLightFullImC4MM

           type(dd_complex) :: result

      result = ddreal(0)

           ampNonresonantLightFullImC4MM = result
       end function ampNonresonantLightFullImC4MM

       function ampNonresonantLightFullImC4MP()
           implicit none
           type(dd_complex) :: ampNonresonantLightFullImC4MP
           type(dd_complex) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           type(dd_complex) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           type(dd_complex) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           type(dd_complex) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           type(dd_complex) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           type(dd_complex) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           type(dd_complex) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           type(dd_complex) ::  t175,t176,t177,t178,t179,t18,t180,t181,t182,t183,t184,t185
           type(dd_complex) ::  t186,t187,t188,t189,t19,t190,t191,t192,t193,t194,t195,t196
           type(dd_complex) ::  t197,t198,t199,t2,t20,t200,t201,t202,t203,t204,t205,t206
           type(dd_complex) ::  t207,t208,t209,t21,t210,t211,t212,t213,t214,t215,t216,t217
           type(dd_complex) ::  t218,t219,t22,t220,t221,t222,t223,t224,t225,t226,t227,t228
           type(dd_complex) ::  t229,t23,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239
           type(dd_complex) ::  t24,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249,t25
           type(dd_complex) ::  t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t26,t260
           type(dd_complex) ::  t261,t262,t263,t264,t265,t266,t267,t268,t269,t27,t270,t271
           type(dd_complex) ::  t272,t273,t274,t275,t276,t277,t278,t279,t28,t280,t281,t282
           type(dd_complex) ::  t283,t284,t285,t286,t287,t288,t289,t29,t290,t291,t292,t293
           type(dd_complex) ::  t294,t295,t296,t297,t298,t299,t3,t30,t300,t301,t302,t303
           type(dd_complex) ::  t304,t305,t306,t307,t308,t309,t31,t310,t311,t312,t313,t314
           type(dd_complex) ::  t315,t316,t317,t318,t319,t32,t320,t321,t322,t323,t324,t325
           type(dd_complex) ::  t326,t327,t328,t329,t33,t330,t331,t332,t333,t334,t335,t336
           type(dd_complex) ::  t337,t338,t339,t34,t340,t341,t342,t343,t344,t345,t346,t347
           type(dd_complex) ::  t348,t349,t35,t350,t351,t352,t353,t354,t355,t356,t357,t358
           type(dd_complex) ::  t359,t36,t360,t361,t362,t363,t364,t365,t366,t367,t368,t369
           type(dd_complex) ::  t37,t370,t371,t372,t373,t374,t375,t376,t377,t378,t379,t38
           type(dd_complex) ::  t380,t381,t382,t383,t384,t385,t386,t387,t388,t389,t39,t390
           type(dd_complex) ::  t391,t392,t393,t394,t395,t396,t397,t398,t399,t4,t40,t400
           type(dd_complex) ::  t401,t402,t403,t404,t405,t406,t407,t408,t409,t41,t410,t411
           type(dd_complex) ::  t412,t413,t414,t415,t416,t417,t418,t419,t42,t420,t421,t422
           type(dd_complex) ::  t423,t424,t425,t426,t427,t428,t429,t43,t430,t431,t432,t433
           type(dd_complex) ::  t434,t435,t436,t437,t438,t439,t44,t440,t441,t442,t443,t444
           type(dd_complex) ::  t445,t446,t447,t448,t449,t45,t450,t451,t452,t453,t454,t455
           type(dd_complex) ::  t456,t457,t458,t459,t46,t460,t461,t462,t463,t464,t465,t466
           type(dd_complex) ::  t467,t468,t469,t47,t470,t471,t472,t473,t474,t475,t476,t477
           type(dd_complex) ::  t478,t479,t48,t480,t481,t482,t483,t484,t485,t486,t487,t488
           type(dd_complex) ::  t489,t49,t490,t491,t492,t493,t494,t495,t496,t497,t498,t5
           type(dd_complex) ::  t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t6,t60
           type(dd_complex) ::  t61,t62,t63,t64,t65,t66,t67,t68,t69,t7,t70,t71
           type(dd_complex) ::  t72,t73,t74,t75,t76,t77,t78,t79,t8,t80,t81,t82
           type(dd_complex) ::  t83,t84,t85,t86,t87,t88,t89,t9,t90,t91,t92,t93
           type(dd_complex) ::  t94,t95,t96,t97,t98,t99

           type(dd_complex) :: result

      t1 = intHLs160000x0112D2eps0()
      t2 = intHLs16s25s26s34s56x1111D2eps0()
      t3 = intHs160000x0112D2eps0()
      t4 = intHs160s26s34s56x1021D2eps1()
      t5 = intHs160s26s34s56x1022D4eps0()
      t6 = intHs16s25s26s34s56x1112D4eps0()
      t7 = intHs16s25s26s34s56x1121D4eps0()
      t8 = intHs16s25s26s34s56x1211D4eps0()
      t9 = intHs16s25s26s34s56x1131D4eps0()
      t10 = ddreal(4)
      t11 = propZ25 * s25
      t12 = -t10 + t11
      t13 = intHs16s25s26s34s56x1221D4eps0()
      t14 = intHs16s25s26s34s56x1311D4eps0()
      t15 = intHs16s25s26s34s56x1311D4eps1()
      t16 = ddreal(2)
      t17 = -t16 + epinv
      t18 = ddreal(3)
      t19 = (gw)**2
      t20 = (gb)**2
      t21 = t18 * s25
      t22 = t21 * t19 * propZ25
      t23 = -t20 * t12
      t24 = t22 + t23
      t25 = ddreal(1) - epinv
      t26 = t20 * t12
      t27 = t26 - t22
      t28 = intHLs16s25s26s34s56x1112D4eps0()
      t29 = t16 + t11
      t30 = t20 * t29
      t31 = t22 + t30
      t32 = intHs16s25s26s34s56x1310D4eps0()
      t33 = s16 + s26 - s34 + s56
      t34 = intHLs16s25s26s34s56x1222D6eps1()
      t35 = s12 + s25
      t36 = t16 * t35
      t37 = t10 * s26
      t38 = t37 + t36 + s16
      t39 = intHLs16s25s26s34s56x1211D2eps1()
      t40 = ddreal(6)
      t41 = ddreal(5)
      t42 = t18 * s26
      t43 = t40 * s12
      t44 = t41 * s25
      t45 = s34 - s56 + t42 + t44 + t43
      t46 = intHLs16s25s26s34s56x1121D2eps1()
      t47 = s25 + s26 + s56
      t48 = intHLs16s25s26s34s56x1112D2eps1()
      t49 = -s13 + s26 + s45 + s46 + s56
      t50 = t16 * t49 + s25
      t51 = intHs160s26s34s56x1031D4eps0()
      t52 = ddreal(1) / s25
      t53 = t20 * (-propZ25 * t16 + ddreal(8) * t52) + t40 * t19 * propZ
     &25
      t54 = intHs16s25s26s34s56x1211D2eps1()
      t55 = intHs16s25s26s34s56x1220D4eps1()
      t56 = intHs16s25s26s34s56x1130D4eps1()
      t57 = t16 * s12
      t58 = t57 + s25
      t59 = intHLs16s25s26s34s56x1122D4eps1()
      t60 = s13 - s45 - s46
      t61 = t16 * s56
      t62 = t61 - t60
      t63 = intHLs16s25s26s34s56x1121D2eps0()
      t64 = intHs16s25s26s34s56x1120D2eps0()
      t65 = intHs16s25s26s34s56x1210D2eps0()
      t66 = intHs16s25s26s34s56x1210D2eps1()
      t67 = intHs16s25s26s34s56x1122D4eps1()
      t68 = t17 * s25
      t69 = t10 * s12
      t70 = t16 * epinv
      t71 = intHLs16s25s26s34s56x1112D2eps0()
      t72 = intHLs16s25s26s34s56x1222D6eps0()
      t73 = intHs16s25s26s34s56x1212D4eps1()
      t74 = t16 * s26
      t75 = intHLs16s25s26s34s56x1211D2eps0()
      t76 = intHLs16s25s26s34s56x1212D4eps1()
      t77 = t16 * t60
      t78 = t18 * s56
      t79 = s25 + s26 - s34 + t69 - t78 + t77
      t80 = epinv * t79
      t81 = intHLs16s25s26s34s56x1212D4eps0()
      t82 = intHs16s25s26s34s56x1222D6eps0()
      t83 = intHs160s26s34s56x1022D4eps1()
      t84 = intHs160000x0112D2eps1()
      t85 = intHs16s25s26s34s56x1121D2eps1()
      t86 = ddreal(1) + epinv
      t87 = t57 * t86
      t88 = intHs16s25s26s34s56x1121D4eps1()
      t89 = intHs16s25s26s34s56x1211D4eps1()
      t90 = intHs16s25s26s34s56x1112D4eps1()
      t91 = -t18 + t70
      t92 = intHLs160000x0111D0eps0()
      t93 = intHLs160000x0111D0eps1()
      t94 = intHLs160000x0112D2eps1()
      t95 = intHLs16s25s26s34s56x1112D4eps1()
      t96 = intHs16s25s26s34s56x1122D4eps0()
      t97 = t74 + s25
      t98 = intHs16s25s26s34s56x1212D4eps0()
      t99 = intHs16s25s26s34s56x1411D6eps0()
      t100 = intHs16s25s26s34s56x1121D2eps0()
      t101 = intHs16s25s26s34s56x1131D4eps1()
      t102 = t68 - t69
      t103 = intHs16s25s26s34s56x1221D4eps1()
      t104 = -t57 + t68
      t105 = intHs16s25s26s34s56x1141D6eps0()
      t106 = intHs16s25s26s34s56x1321D6eps0()
      t107 = t57 + t21
      t108 = intHs16s25s26s34s56x1231D6eps0()
      t109 = t69 + t21
      t110 = intHLs16s25s26s34s56x1132D6eps0()
      t111 = intHLs16s25s26s34s56x1312D6eps0()
      t112 = s25 + s16
      t113 = t16 * (s12 + s26)
      t114 = t112 + t113
      t115 = intHs16s25s26s34s56x1312D6eps0()
      t116 = t21 + t113
      t117 = intHs16s25s26s34s56x1113D4eps1()
      t118 = intHs16s25s26s34s56x1114D6eps0()
      t119 = intHs16s25s26s34s56x1213D6eps0()
      t120 = s12 + s25 + s26
      t121 = t120 * t16 + s34
      t122 = intHLs16s25s26s34s56x1113D4eps0()
      t123 = -t78 + t77 - s26
      t124 = intHLs16s25s26s34s56x1122D4eps0()
      t125 = intHLs16s25s26s34s56x1114D6eps0()
      t126 = s25 + s26 - s56
      t127 = intHLs16s25s26s34s56x1123D6eps0()
      t128 = t16 * s25
      t129 = t42 + t128 - s56
      t130 = intHLs16s25s26s34s56x1213D6eps0()
      t131 = t42 + t36 - s56
      t132 = intHs16s25s26s34s56x1132D6eps0()
      t133 = s26 + s34 - s56 + t69 + t128
      t134 = intHs16s25s26s34s56x1123D6eps0()
      t135 = t16 * (s12 + s34) + t126
      t136 = intHs16s25s26s34s56x1310D4eps1()
      t137 = intHs16s25s26s34s56x1110D2eps0()
      t138 = intHs160s26s34s56x1020D2eps0()
      t139 = intHLs16s25s26s34s56x1113D4eps1()
      t140 = intHs160s26s34s56x1031D4eps1()
      t141 = intHLs16s25s26s34s56x1123D6eps1()
      t142 = intHLs16s25s26s34s56x1213D6eps1()
      t143 = intHs16s25s26s34s56x1132D6eps1()
      t144 = intHs16s25s26s34s56x1123D6eps1()
      t145 = intHs16s25s26s34s56x1110D2eps1()
      t146 = intHs16s25s26s34s56x1220D4eps0()
      t147 = intHs16s25s26s34s56x1130D4eps0()
      t148 = intHs16s25s26s34s56x1120D2eps1()
      t149 = t25 * s25
      t150 = intHs160s26s34s56x1020D2eps1()
      t151 = intHs16s25s26s34s56x1321D6eps1()
      t152 = intHs16s25s26s34s56x1231D6eps1()
      t153 = intHLs16s25s26s34s56x1132D6eps1()
      t154 = intHLs16s25s26s34s56x1312D6eps1()
      t155 = intHs16s25s26s34s56x1312D6eps1()
      t156 = intHs16s25s26s34s56x1213D6eps1()
      t157 = ddreal(1) / t33
      t158 = t17 * t145
      t159 = t150 * t25
      t160 = t25 * t34
      t161 = t160 - t72
      t162 = t48 * epinv
      t163 = epinv * t39
      t164 = t46 * t47
      t165 = t164 * epinv
      t166 = t161 * t38
      t167 = (t162 + t71) * t50
      t168 = t25 * t66
      t169 = (-t159 + t138) * t52
      t170 = (t148 * (t57 + t149) + t137 + t158) * t52 * t24
      t171 = t67 * (t70 * s26 - s26 - s34 + s56 + t68 - t69)
      t172 = t25 * intHs16s25s26s34s56x1222D6eps1()
      t173 = t88 + t89
      t174 = t7 + t8
      t79 = t79 * t81
      t175 = t47 * t63
      t176 = t91 * t95
      t177 = t84 * epinv
      t178 = t85 * (t87 + s25)
      t179 = t83 * t25
      t180 = t1 * t29
      t181 = t90 * t91
      t182 = t25 * intHs16s25s26s34s56x1114D6eps1()
      t183 = t25 * intHLs16s25s26s34s56x1114D6eps1()
      t184 = t183 - t125
      t185 = t25 * intHs16s25s26s34s56x1141D6eps1()
      t186 = -t185 + t105
      t187 = t9 + t13 + t14
      t188 = intHs16s25s26s34s56x1411D6eps1() * t25
      t189 = (-t182 + t118) * s34
      t190 = t186 * t58
      t191 = t187 * s25
      t192 = t184 * t31
      t193 = t25 * t55
      t194 = t25 * t56
      t195 = t147 - t194
      t196 = t97 * t153
      t197 = t114 * t154
      t198 = t107 * t151
      t199 = t116 * t155
      t200 = t135 * t144
      t201 = t195 * t58
      t202 = (t146 - t193) * t35
      t203 = t15 * t17
      t204 = t52 * t31
      t205 = t25 * t52
      t206 = s34 * t117
      t207 = t108 * t109
      t208 = t115 * t116
      t107 = t106 * t107
      t209 = t134 * t135
      t210 = t101 * t102
      t211 = t140 * t25
      t212 = t100 * s12
      t213 = t25 * t136
      t214 = t122 * t123
      t215 = -t124 * t62
      t216 = t97 * t110
      t114 = t111 * t114
      t217 = t2 * t29
      t129 = t24 * (t52 * (t103 * t104 + t210 - t211 + t212) - t213 * t1
     &57) + t24 * (t157 * (t52 * (t202 + t201) + t32) + t203) - t26 * t1
     &87 + t204 * (-t130 * t131 + t25 * (t129 * t141 + t131 * t142 + t19
     &6 + t197) - t28) + t205 * (-t109 * t152 - t121 * t156 - t133 * t14
     &3 - t198 - t199 - t200) * t27 + (t119 * t121 + t132 * t133 + t107 
     &+ t206 + t207 + t208 + t209) * t52 * t27 + t204 * (-t127 * t129 - 
     &t139 * (-epinv * t123 + s25 + s26 - s56) + t214 + t215 - t216 - t1
     &14) + t52 * (t12 * t6 + t217) * t20
      t131 = t19 * propZ25
      t187 = intHL0s25s26s34s56x1130D4eps0()
      t218 = -s13 + s25 + s26 + s56
      t219 = s45 + s46
      t220 = t16 * t219
      t221 = t218 + t220
      t222 = s13 - s45 - s56
      t223 = -t16 * t222 + s25 + s46
      t224 = intHLs160000x0211D2eps0()
      t225 = -s16 + s25 + s34
      t226 = t16 * (s12 + s13 - s56)
      t227 = t18 * s46
      t228 = t10 * s45
      t229 = t225 + t228 + t227 - t226
      t230 = t10 * t219
      t226 = t225 + t230 - t226
      t231 = intHs160000x0121D2eps0()
      t232 = s25 - s26 - s16 + s34
      t233 = s13 - s56
      t234 = t18 * t233
      t235 = -t57 + t232 + t230 - t234
      t236 = intHLs16s25s26s34s56x1211D4eps0()
      t237 = s25 - s26
      t238 = s16 - s34
      t239 = t16 * t238
      t240 = t40 * s45
      t241 = t41 * s46
      t242 = t237 - t69 - t234 - t239 + t241 + t240
      t243 = intHs16s25s26s34s56x1111D2eps0()
      t244 = t16 * (s12 - s45 - s46)
      t245 = s26 + s16 - s34 + s13 - s56 + t244
      t246 = intHL0s25s260s56x1012D2eps0()
      t247 = t47 - t77
      t248 = intHL0s25s260s56x1013D4eps0()
      t249 = intHL0s25s260s56x1021D2eps0()
      t250 = intHs160000x0111D0eps0()
      t251 = intHL0s25s26s34s56x1210D2eps1()
      t252 = intHL0s25s26s34s56x1310D4eps1()
      t253 = t131 * s25
      t254 = t253 + t20
      t255 = s25 + s16 - s34
      t256 = t255 + t113
      t257 = intHs160s26s34s56x1013D4eps1()
      t258 = s26 + s16 - s34 - s56
      t259 = t16 * (s12 + s13 - s45 - s46)
      t260 = t258 + t259
      t261 = t16 * t33
      t262 = t261 + s25
      t263 = s45 + s46 + s56
      t264 = t16 * t263
      t265 = t264 + s25
      t266 = t57 * s13
      t267 = t33 - t128
      t268 = s13 * (-t21 + t261) - t265 * t267
      t269 = s25 + s56
      t270 = s26 + s56
      t271 = t16 * t270
      t272 = t16 * s13
      t273 = t265 * (t21 + t271) - t272 * (t16 * t269 + s26)
      t274 = (t10 * t33 + t57) * s13 - t265 * (t261 - s25)
      t275 = intHL0s25s260s56x1022D4eps0()
      t276 = t270 + t128
      t277 = s13 * t276
      t278 = t270 * t263
      t264 = (-t264 - s26 - s25) * s25 + t277 - t278
      t279 = s12 + s13
      t280 = t74 + t255
      t281 = t279 * s25
      t282 = propW16 * t256
      t283 = t282 * t47 * t33
      t284 = t10 * t279
      t285 = s25 + s26
      t286 = s12 + s13 + s16
      t287 = s12 + s16 + s26 - s34
      t288 = t286 * s25
      t289 = (-s13 + s25 + s26) * s34
      t290 = t219 * t225
      t291 = (s12 + s16 + s25 + s26) * s26
      t292 = s13 * s16
      t293 = (s56)**2
      t294 = s12 + s13 - s16
      t295 = s25 - s12
      t296 = t18 * s13
      t297 = -t295 + t296
      t298 = -t272 + s25
      t299 = s16 - s34 - s13
      t300 = -t16 * t299 - t295 - t42
      t301 = t10 * t60 - t16 * t255 + s12 - t42
      t302 = s12 - s26
      t303 = s13 - s45
      t304 = t279 * s13
      t305 = (-s34 + t286) * s25 + (s25 - s12 - s13) * s26 + (-t272 - s2
     &5 - s12 + s26 + s45) * s45 + (-t16 * t303 + s46 - t302) * s46 + (-
     &t302 - t77 + s56) * s56 + t304
      t306 = t294 * s25
      t307 = (-t16 * (s25 + s16 - s34 - s13 + s45) - s12 - t42 - s46) * 
     &s46
      t308 = t272 * s16
      t309 = t18 * t293
      t310 = s26 + s34
      t311 = s25 * t219
      t312 = (s25)**2
      t313 = t10 * t311
      t314 = t16 * (s34 * t263 + t312) + s13 * (-t21 + t239) - s16 * t26
     &5 + t310 * s25 + t313 + t21 * s56
      t315 = s13 * (t61 + s25)
      t316 = t47 * t265
      t317 = -t316 + t315
      t318 = intHL0s25s260s56x1031D4eps1()
      t319 = t16 * (s13 - s45 - s46 - s56)
      t320 = -t319 + s25
      t321 = s26 - s56
      t322 = s46 * t321
      t323 = intHLs16s25s26s34s56x1311D4eps1()
      t324 = t16 * s45
      t325 = t324 + s25 + s26 - s13 + s46 + s56
      t326 = -s16 + s25
      t327 = t18 * s34
      t328 = -t326 - t327
      t329 = t47 + t220
      t330 = s26 * s34
      t331 = s13 * t328 - s16 * t329 + s25 * t47 + t10 * s34 * t219 + t1
     &6 * (s34 + s45 + s46) * s25 + t330 + t327 * s56
      t332 = intHL0s25s26s34s56x1220D4eps0()
      t333 = t16 * s34
      t334 = t47 + t333
      t335 = s25 + s26 + s34 + s56
      t336 = s13 * t334 - t329 * t335
      t337 = intHL0s25s26s34s56x1120D2eps1()
      t338 = t70 * s34
      t339 = epinv * s34
      t340 = t339 - s25 - s26 - s56
      t341 = -t16 * (s12 + s26 + s13 - s45 - s46 - s56) - t238
      t342 = ddreal(7)
      t343 = t40 * t219
      t344 = t342 * s56
      t345 = -s26 - s45 - s56
      t346 = s26 * s46
      t347 = s46 * s56
      t348 = t41 * s13
      t349 = t348 * s25
      t350 = -t16 * (s13 * t270 - s26 * s45 + s56 * t345) + t18 * (t347 
     &+ t312) + s25 * (t344 + t343 + s26) + t346 - t349
      t351 = intHL0s25s260s56x1021D2eps1()
      t352 = t25 * s26
      t353 = epinv * s56
      t354 = -s34 * t298
      t355 = (-t300 + s45) * s45
      t356 = -s26 * t297 - s56 * t301 + t304 - t306 - t307 - t308 + t309
     & + t354 + t355
      t357 = intHLs16s25s26s34s56x1231D6eps1()
      t358 = s13 + s16
      t359 = t57 + t358
      t360 = -s13 + s16
      t361 = s25 - s13
      t362 = t16 * t285
      t363 = s25 + s45 + s46
      t302 = t16 * t302
      t364 = t10 * t363
      t365 = (t57 + t360 + s26) * s26
      t366 = (t359 - s25) * s25
      t367 = s12 * s13
      t368 = t16 * (t219 * (t362 + s34) + t367) + t309 + (t16 * t361 + s
     &26) * s34 + (t18 * (s34 - s13) - s16 - t302 + t364) * s56 + t292 -
     & t365 - t366
      t369 = t18 * s16
      t370 = t41 * s34
      t371 = t342 * s25
      t372 = (t57 + s25 + s16 + s13 + s26) * s26
      t310 = t10 * (t310 + t128) * t219 - t16 * (-s16 * t60 - t312 - t36
     &7) + t309 - (t284 + t369) * s25 + (t10 * t361 + s26) * s34 + (-t18
     & * t358 + t230 - t302 + t370 + t371) * s56 - t372
      t373 = s16 + s26 + s56
      t374 = t373 - t327
      t375 = t373 - t333
      t220 = t258 - t220
      t376 = s13 * t374 + t220 * t375 + t57 * t375
      t377 = t25 * s34
      t378 = -t377 + s26 + s16 + s56
      t379 = t33 + t338
      t380 = intHs160000x0211D2eps1()
      t381 = t18 * epinv
      t382 = t381 - t10
      t383 = intHLs16s25s26s34s56x1321D6eps1()
      t384 = s25 + s16 - s13
      t385 = t16 * t295
      t386 = t10 * s34
      t387 = -t16 * (t219 * (t285 + t333) + t367) + (t57 + s16) * s25 + 
     &(t57 + t384 + s26) * s26 + (t296 - t385 + s16 - s34) * s34 + (-s25
     & + s16 + s13 + t244 - t386 - s56) * s56 - t292
      t388 = t10 * s13
      t389 = s26 + s16
      t390 = s26 + s45 + s46
      t391 = t360 * s25
      t392 = (t74 - t388 + t21) * s34
      t393 = -s12 - s16
      t394 = t393 * s13
      t384 = -t16 * (t219 * (t389 - t327) + t394) - (t112 - t296 + s26) 
     &* s26 - (t16 * t390 + s56 + t384 - t386) * s56 - t391 + t392
      t395 = s26 - s13
      t396 = t16 * t395
      t397 = t16 * (s26 - s13 + s45)
      t398 = s25 - s12 + s16 - s34
      t399 = t361 * s34
      t400 = s16 * s25
      t294 = t16 * t293 + s13 * t294 + (-t272 + t295) * s26 + (t396 - s1
     &2 + s16 - s34 + s45) * s45 + (t397 + t255 + s46) * s46 + (-t18 * t
     &60 + t398 + t74) * s56 - t399 + t400
      t295 = t16 * s16
      t401 = (s34)**2
      t303 = -t16 * (s16 * t279 - (s12 + s13 + s16 + s26) * s34 - t293) 
     &- (t18 * t279 + s26 + t295) * s26 - (-t397 - t398 - s46) * s46 + (
     &-t10 * t303 - t18 * (s12 - s46) + t285) * s56 - t281 + t304 + t355
     & - t401 - (s16)**2
      t355 = -t326 - t333
      t328 = t16 * t355 * t219 - (t57 + s13 + s16) * s16 + s25 * t359 + 
     &s26 * t326 - (t16 * (s25 - s12 - s16) - s26 - t296 + s34) * s34 + 
     &s56 * t328
      t359 = -t393
      t398 = t16 * t359
      t402 = t10 * s25
      t399 = t16 * (t292 + t399) - (s26 + t398 + t296 - s25) * s25 + (t4
     &02 - t239) * s45 - (-t402 + t239 + s26) * s46 + (t21 - t239 + s46)
     & * s56
      t403 = t348 - t128 - s16
      t404 = t18 * t360
      t390 = t10 * t390
      t360 = t128 * t360
      t405 = -t10 * t394 - t309 + (t403 - s26) * s26 + (t21 - t348 + s26
     &) * s34 + (-t390 - t404 + t370 - t128) * s56 + t219 * (s34 * t40 -
     & t10 * t389) - t360
      t406 = -s13 * t342 + s26 + t402
      t407 = t342 * s34
      t237 = t16 * (t367 + t312) - t309 - (t369 + s13) * s25 - (t112 - t
     &348 + s26) * s26 + t406 * s34 + (-s16 * t41 + s25 + t296 - t390 + 
     &t407) * s56 + t219 * (-s16 * t40 + t10 * t237 + ddreal(8) * s34) +
     & t348 * s16
      t408 = intHs16s25s26s34s56x1113D4eps0()
      t409 = intHLs160000x0211D2eps1()
      t410 = intHs160000x0121D2eps1()
      t411 = intHL0s25s260s56x1013D4eps1()
      t412 = intHL0s25s260s56x1010D0eps0()
      t413 = intHL0s25s260s56x1011D2eps0()
      t414 = intHL0s25s260s56x1020D2eps0()
      t415 = intHL0s25s26s34s56x1110D2eps0()
      t416 = intHL0s25s26s34s56x1130D4eps1()
      t417 = intHLs160000x0111D2eps0()
      t418 = intHLs16s25s26s34s56x1411D6eps0()
      t419 = intHs160s26s34s56x1011D2eps0()
      t420 = t16 * (s12 + s13 - s45) + s26 + s16 - s34 - s46 - s56
      t421 = intHLs16s25s26s34s56x1311D4eps0()
      t422 = s26 - s16 + s34 - t234 + t343 + t385
      t423 = intHLs16s25s26s34s56x1121D4eps0()
      t424 = s46 + s56
      t240 = t41 * t424 + s26 + t21 + t240 - t388
      t425 = s16 + s26 - s34
      t343 = s56 * t41 + t343 + t385 - t388 - t425
      t241 = -t222 * t40 + t21 + t241
      t385 = s26 - s34
      t426 = (s46)**2
      t286 = -t16 * (s13 * t286 + (s25 - s16 + s34 - s13) * s45 - t288 +
     & t289 - t291 - t426) + (t16 * (s16 + s45) + t385) * s46 + (t16 * (
     &s12 + s26 + s16 - s34 + s13) + t227) * s56
      t427 = t97 * s34
      t395 = -t395 - t128
      t428 = (s13)**2
      t429 = t285 * s34
      t430 = (s12 + s13 + s16 + s25 + s26) * s26
      t431 = (t272 - s25 + s12 - s45 + s46 - s56) * s56
      t324 = (-t324 + s25) * s46
      t432 = (t69 + t348 + t295) * s25
      t433 = intHLs16s25s26s34s56x1321D6eps0()
      t434 = intHLs16s25s26s34s56x1131D4eps0()
      t279 = t16 * t279
      t435 = s12 - s45
      t436 = t16 * t435
      t437 = (t279 + s16) * s25
      t438 = intHLs16s25s26s34s56x1221D4eps0()
      t439 = -t42 - t128 + s13
      t440 = s25 + s26 + s34
      t441 = -s25 + s16 - s34 + s13 - s46 + t436
      t442 = t440 * s46
      t443 = t16 * (s45 * t440 - t367)
      t444 = (s12 - s13 + s16 + s25 + s26) * s26
      t396 = (-t396 - s25) * s34 + t16 * ((t16 * t389 + s12 - s13 + s25 
     &- s34 + s45 + s46 + s56) * s56 + t219 * (t385 + t295) - t367 + t40
     &0 + t444) - t388 * s16
      t400 = t16 * (s26 + s16 - s13 + s45 + s46)
      t333 = t219 * (t389 - t333)
      t401 = t16 * (-t292 + t293 - t401 + t333) + (t358 + s26) * s25 + (
     &t284 + t42 + t295) * s34 + (-t327 + t400 + s25) * s56
      t445 = s25 - s34
      t446 = t298 * s26
      t447 = t219 * t425
      t358 = t16 * (-t292 + t293 + t447) + s25 * t358 + (t389 + t279 - s
     &34) * s34 + (t445 + t400) * s56 + t446
      t389 = t296 - t128 - s26
      t299 = (t16 * (s12 + s45 + s46) + t18 * t299 + s25 + t37) * s56
      t361 = (t18 * t361 + s16 + s26 + t57) * s26
      t232 = t219 * t232
      t400 = t16 * (t367 - t288 + t232)
      t448 = t18 * (t292 - t293)
      t449 = intHs16s25s26s34s56x1211D2eps0()
      t450 = intHLs160000x0121D2eps0()
      t451 = t282 * s16
      t452 = t16 * (t451 - s25)
      t453 = t451 * t40
      t454 = t10 * t222 - t227 + t452
      t227 = t29 * (-t10 * (s45 + s56) - t128 - t227 + t388) + t453
      t279 = (t369 + t279 - s25) * s25
      t444 = t16 * ((-s25 + s12 + s26 + t239) * s56 - t367 + t444)
      t455 = t10 * (t292 + t290)
      t456 = t388 - t21
      t457 = s26 + s16 - s34 - s13 + s45 + s46
      t458 = -s34 * t389
      t459 = intHLs16s25s26s34s56x1221D4eps1()
      t460 = t16 * t457 + s25
      t447 = t293 + t394 + t447
      t461 = t16 * t298
      t435 = t435 * s25
      t462 = (-t345 - s13) * s13
      t258 = t258 * s25
      t346 = -t16 * (s45 * s46 + t426 + t435 + t462) - t18 * (s13 * s25 
     &+ t347) - t346 - t258
      t347 = s25 * s26
      t258 = -t16 * (-s46 * t345 + t426 + t435 + t462) - t258 - t21 * s1
     &3
      t345 = intHLs16s25s26s34s56x1141D6eps0()
      t435 = (s26)**2
      t462 = s12 * t47
      t463 = t363 * s16
      t464 = intHLs16s25s26s34s56x1131D4eps1()
      t465 = intHLs160000x0121D2eps1()
      t466 = intHs160s26s34s56x1012D2eps1()
      t467 = intHs16s25s26s34s56x1112D2eps1()
      t468 = intHL0s25s26s34s56x1220D4eps1()
      t469 = intHL0s25s26s34s56x1110D2eps1()
      t470 = intHL0s25s260s56x1010D0eps1()
      t471 = intHL0s25s260s56x1020D2eps1()
      t472 = intHL0s25s26s34s56x1120D2eps0()
      t473 = intHL0s25s26s34s56x1210D2eps0()
      t474 = intHL0s25s26s34s56x1310D4eps0()
      t475 = intHLs160000x0122D4eps0()
      t476 = intHs160s26s34s56x1011D2eps1()
      t477 = intHs160s26s34s56x1013D4eps0()
      t478 = intHL0s25s260s56x1011D2eps1()
      t479 = intHL0s25s260s56x1022D4eps1()
      t480 = intHLs16s25s26s34s56x1231D6eps0()
      t481 = t18 * t238
      t482 = intHs160000x0211D2eps0()
      t483 = intHL0s25s260s56x1031D4eps0()
      t484 = intHLs16s25s26s34s56x1211D4eps1()
      t485 = t10 * t233
      t486 = ddreal(8) * s45
      t487 = t342 * s46
      t488 = ddreal(8) * t219
      t489 = intHLs16s25s26s34s56x1121D4eps1()
      t490 = ddreal(1) / t47
      t491 = t25 * intHL0s25s260s56x1012D2eps1()
      t492 = s34 * t75
      t493 = t89 * t382
      t494 = intHs160000x0111D0eps1() * epinv
      t495 = t24 * t52
      t496 = ddreal(1) / (ecossin)**2 * propW34 * im
      t39 = s34 * t39
      t497 = s34 * intHs16s25s26s34s56x1112D2eps0()
      t498 = t226 * t3 + t231 * t235 + t243 * (s13 + s26 - s56) + t96 * 
     &t401 + t497 * t320
      t297 = t24 * t498 + t27 * (t100 * t396 + t237 * t82 + t343 * t7 + 
     &t422 * t6 + t148 * (-t309 + epinv * t305 + s26 * t297 + s34 * t298
     & + (t300 - s45) * s45 + s56 * t301 - t304 + t306 + t307 + t308) * 
     &t157) + t31 * (-t124 * (t16 * (s45 * t395 + t426 - t428 - t429 + t
     &430 + t431) + t432 - t324) - t2 * (t57 + s25 + s26 + s16 - s34 - s
     &13 - s46 + s56) + t223 * t92 - t224 * t229 + t236 * t242 + t28 * t
     &241 + t63 * (-t16 * ((s25 + s26 - s16) * s45 - t288 - t291 + t292 
     &+ t293) + (-t285 + t295) * s46 + (-t16 * (s25 - s12 - s16 - s13 + 
     &s45) - s46) * s56 - t427) - t81 * t286 + t423 * t240 - t438 * (s34
     & * t439 + (t441 - s56) * s56 + t372 + t437 - t442 - t443) + (t337 
     &* (s13 * (t47 - t338) + t329 * t340) - t25 * t336 * t468) * t490 -
     & t39 * (epinv * t223 + s25 - t319))
      t298 = t496 * t52
      t300 = intHs160s26s34s56x1012D2eps0() + intHs160s26s34s56x1021D2ep
     &s0()
      t301 = t157 * t25
      t304 = -t260 * t157
      t306 = s34 * t251
      t307 = (t472 + t473) * s34
      t43 = t496 * (t52 * (t24 * ((t86 * t320 * t467 + t304 * t300) * s3
     &4 + t157 * (t145 * t420 + t260 * t476) * t17) + t27 * (t221 * t250
     & + t88 * (epinv * t343 + s16 + s26 - s34 - t21 - t344 + t348 - t48
     &8 + t57) + t90 * (epinv * t422 + s16 - s34 - t21 + t302 + t485 - t
     &488) + t301 * (t150 * t303 * t157 - t305 * t66 - t376 * t83)) + t3
     &1 * (t484 * (epinv * t242 - s25 + t43 + t481 + t485 - t486 - t487 
     &+ t74) + t95 * (epinv * t241 + ddreal(8) * t222 - t402 - t487) + t
     &489 * (epinv * t240 - t342 * t424 - s26 + t348 - t402 - t486))) + 
     &t31 * (t52 * (t17 * (-t469 - t478) + t25 * (-t306 - t470 + t471) +
     & t307) + t246 + t249) * t490 * t247)
      t66 = epinv * t465
      t145 = t52 * (t103 * (t10 * (t292 - t293 + t232) + t266 + epinv * 
     &(-t361 - t299 + t458 + t448 + t400) - (t57 + t369 + s13 - s25) * s
     &25 + s26 * t456 - t456 * s34 + (-t10 * t457 + s25) * s56) + t73 * 
     &(epinv * (t16 * (t367 - t435 + t290) - (t398 + s13) * s25 - (t16 *
     & (s12 - s13) + s16 + t21) * s26 + (-t481 - t113 + s25) * s56 + t45
     &8 + t296 * s16) - t279 + t392 + t455 - t444)) - t54 * (epinv * t22
     &1 - s13 - s16 - s26 + s34 + s56 - t244)
      t222 = t52 * (t66 + t450)
      t232 = t204 * (t310 * t72 - t46 * (epinv * (-t16 * (s13 * (s25 - s
     &16 + s56) - s25 * s45 + (s25 + s16 - s34 - s45) * s26 + (-s25 + s1
     &6 - s45) * s56 - t293 + t435 + t462 + t463) + s25 * s34 + s46 * t4
     &7) - t315 + t316) - t48 * (epinv * t258 + s25 * t221) + t59 * (t10
     & * (t219 * (t128 + s26) + t293 + t312) - t349 + epinv * (-t16 * t3
     &95 * s45 - t16 * (t431 + t430 - t429 + t426 - t428) + t324 - t432)
     & + (-t272 + t21) * s26 + (ddreal(9) * s25 + t10 * (s26 - s13 + s45
     & + s46)) * s56) - t71 * t258 + t76 * (-epinv * t286 - t279 + t392 
     &- t444 + t455) + t459 * (t10 * t440 * t219 + t266 + t309 + epinv *
     & (-s34 * t439 + (-t441 + s56) * s56 - t372 - t437 + t442 + t443) +
     & (-t388 + t21 + s26) * s34 + (t364 - t296 + t370 - t302 - s16) * s
     &56 + t292 - t365 - t366))
      t145 = t496 * (t24 * (t221 * t449 + t52 * (t13 * (s34 * t389 + t29
     &9 + t361 - t400 - t448) + t67 * (t10 * (t394 + t333) + t309 + t360
     & + epinv * t401 + (-t403 + s26) * s26 - t406 * s34 + (t390 - t407 
     &+ t404 + t128) * s56) - t85 * (epinv * t396 + s56 * t460 + t16 * t
     &447 + t354 + t391 + t446) + t98 * (t16 * (s12 * t218 + s26 * t270 
     &- s34 * t363 - t311 + t463) + t18 * (s56 * t238 + t347) + s13 * (-
     &t74 - t481 + s25) - s25 * s56 + s26 * t238) - t377 * t304 * t466))
     & + t27 * t145 + t52 * (t1 * (t20 * t227 + t22 * t454) + t94 * (t20
     & * (epinv * t227 - t453) + t253 * (t381 * t454 - t453))) + t222 * 
     &(t20 * (t29 * (-t18 * t424 - s26 - t128 - t228 + t296) + t453) + t
     &22 * (t18 * (s13 - s46 - s56) - s26 - t228 + t452)) + t232)
      t218 = epinv * t380
      t227 = t25 * t252
      t228 = t151 * t25
      t232 = t101 * t24
      t240 = t25 * t257
      t60 = t52 * (t24 * (-t25 * (t136 * t305 + t294 * t55 + t356 * t56)
     & + (-t240 + t477) * t260 * s34) + epinv2 * (t20 * (t16 * (t219 * (
     &t238 + t42 + t128) - t292) + t18 * (t283 + t293) - (-propZ25 * (s5
     &6 * t287 + t288 - t289 - t290 + t291 - t292) + s16 + t284) * s25 -
     & (s13 * t40 + s26 + t112 + t69) * s26 + (t272 + t285) * s34 - (t40
     & * t60 - t21 + t238 + t69 - t74) * s56) - t22 * (s13 * t238 + (t27
     &2 + s12) * s26 + (-s25 + s12 - s26 + t77 - s56) * s56 - t219 * t28
     &0 + t281 - t283)) * t490) * t157
      t60 = t496 * (t31 * (t52 * (-t111 * t328 + t127 * t350 + t130 * t3
     &99 + t139 * (-t18 * s25 * t233 + epinv * t346 + t16 * t312 + t313 
     &- t322 + t347) - t187 * t221 - t387 * t433 - t434 * (-t362 * s45 -
     & s46 * t285 + (-s25 + s16 + s13 - s46 + t436 - s56) * s56 + t372 -
     & t427 + t437) + t490 * (-t25 * t264 * t479 + (t227 - t474) * t247 
     &* s34)) - t248 * t247 * t490) + t52 * (t27 * (t115 * t314 + t119 *
     & t331 + t132 * t405 + t134 * t384 - t228 * t268 + t245 * t51) + t2
     &32 * (t16 * t460 * s56 + epinv * t358 + s26 * t461 - t461 * s34 + 
     &t10 * t447 + t360)) + t60)
      t112 = t496 * (t27 * (t14 * t320 + t15 * (epinv * t320 - s16 + s25
     & - s26 + s34 + t230 - t234 - t57) + t52 * (t106 * t268 + t108 * t2
     &74)) + t204 * (t110 * t273 + t122 * t346 - t320 * t483 + t368 * t4
     &80 + t464 * (-epinv * (t16 * (-s45 * t47 - t330 + t462) + (s16 - s
     &34 - s46) * s25 + (s25 + s16 - s46) * s26 + (-s25 + s16 - s46) * s
     &56 + t277 - t293 + t435) - t16 * t317)) + t495 * (t320 * t482 + t3
     &58 * t9))
      t60 = t112 + t60 + t496 * (t25 * (t31 * (t52 * (-t141 * t350 - t14
     &2 * t399 - t153 * t273 + t154 * t328 + t221 * t416 + t318 * t320 -
     & t357 * t368 + t383 * t387) + t247 * t411 * t490) - (t140 * t245 +
     & t143 * t405 + t144 * t384 + t152 * t274 + t155 * t314 + t156 * t3
     &31) * t52 * t27) + t52 * (t31 * ((t323 * (epinv * t325 - s16 + s25
     & - s26 + s34 + t230 - t234 - t57) + t325 * t421) * s34 + t275 * t2
     &64 * t490) + s34 * (t117 * (epinv * t341 + s26 + t128 + t230 - t23
     &4) + t341 * t408) * t27) + t495 * (t157 * (t146 * t294 + t147 * t3
     &56 + t305 * t32) + t218 * t320))
      t112 = t25 * intHLs16s25s26s34s56x1141D6eps1()
      t219 = t112 - t345
      t230 = t25 * intHLs16s25s26s34s56x1411D6eps1()
      t233 = (t230 - t418) * s34
      t234 = t25 * intHLs160000x0122D4eps1()
      t241 = intHLs160000x0111D2eps1() * t17
      t242 = t496 * (t52 * (t27 * ((s13 * t262 - t265 * t33 + t266) * (-
     &t185 + t105) + (-t182 + t118) * t221 * s34) + t31 * (t317 * t219 +
     & (-t183 + t125) * (t16 * s25 * (t263 - s13) + t312 - t322) + t233 
     &* t245)) + t282 * t52 * ((t234 - t475) * s16 - t417 - t241) * t254
     & + t27 * (-t188 + t99) * t320)
      t244 = -t25 * intHLs160000x0113D4eps1() + intHLs160000x0113D4eps0(
     &)
      t245 = ddreal(1) / ddreal(3)
      t258 = ddreal(2) / ddreal(3)
      t263 = -s16 + s25 - s26 + s34 - s56
      t264 = s25 + s34
      t265 = epinv * t33
      t266 = epinv * s25
      t268 = t26 * t14
      t232 = t24 * (-t15 * t86 + t157 * (t25 * (s34 * t257 * t52 + t136)
     & - t32) - t115 * t225 * t52 - t52 * (t106 + t9) * t263) + t31 * ((
     &-t52 * t94 + t464) * epinv + t490 * (t25 * (s34 * t252 * t52 + t41
     &1) - t248) + t52 * (s34 * t111 + t110 * t47 + t127 * t276 + t130 *
     & t225 + t139 * (t339 + s25)) + t434) + t52 * (-t27 * t33 * t108 - 
     &t180 * t20) + t52 * ((-t232 * t263 - t27 * t380) * epinv + (t24 * 
     &(-t117 * t86 - t157 * t477 - t119 - t408) + t31 * (-t474 * t490 + 
     &t122 + t323)) * s34 + t20 * (-t12 * (t147 + t482) - t29 * t483) + 
     &t25 * (t24 * (s34 * t156 + t151 * t263 + t155 * t225) + t27 * (t15
     &2 * t33 + t56) + t31 * (-s34 * t154 - t141 * t276 - t142 * t225 - 
     &t153 * t47 + t318))) + t268
      t252 = t17 * t476
      t257 = t146 * t263
      t272 = t148 * (t266 + s26 + s16 - s34 + s56)
      t273 = t300 * s34
      t274 = t470 - t471
      t275 = t275 * t276
      t277 = t351 * (t149 + s26 + s56)
      t279 = (t469 + t478) * t17
      t281 = (s34 * t5 + t252 + t257) * t24
      t282 = (t273 + t272) * t27
      t283 = s34 * t467
      t284 = epinv * t85
      t286 = t25 * t88
      t288 = epinv * t410
      t267 = t495 * (-t103 * (epinv * t267 + s16 + s26 - s34 + s56) - t1
     &3 * t267 + t157 * t419 + t225 * t98 - t96 * t355 + t67 * (-epinv *
     & t355 + s34) + t493 - t286)
      t289 = t95 * t382
      t290 = t276 * t479
      t291 = t490 * (-t246 - t249) + t52 * (s34 * t71 + t326 * t63 - t28
     &9 - t423) + t205 * (-t290 * t490 + t335 * t34 + t489) + t52 * (-s3
     &4 * t337 * t490 + t465) * epinv
      t292 = (s34 * t83 + t263 * t55) * t25
      t293 = (t64 + t65 - t205 * (s34 * t466 + t150)) * t27
      t23 = t23 * t52
      t294 = t52 * t27
      t239 = t157 * (t495 * (t158 - t292) + t293) + t31 * t291 + t31 * (
     &t52 * (-s34 * t81 - t124 * t264 + t33 * t438 - t335 * t72 + t459 *
     & (t265 - s34) + t490 * (-t412 + t413 + t414 + t415) - t59 * (epinv
     & * t264 + s25 + s26 + s56) - t76 * (s34 * t86 - s16 + s25)) + t491
     & * t490) + t52 * (t27 * (t375 * t82 + t288) + t30 * t450) + t52 * 
     &(t157 * (-epinv2 * (t20 * ((-propZ25 * t225 + t10) * s25 + t270 * 
     &t40 + t239) + t22 * (t255 + t271)) * t490 + t281 + t282) + t31 * (
     &(s34 * t48 + t326 * t46) * epinv + t490 * ((t332 - t472 - t473) * 
     &s34 + t25 * t274 + t275 + t277 + t279) + t39 * t86) + t86 * (t24 *
     & t73 * t225 + t283 * t27) + t26 * t231 + t284 * t326 * t27) + t267
     & + t23 * (t243 + t7) + t294 * (t100 * t326 - t172 * t375 + t497)
      t255 = t24 * t137
      t264 = t168 * t27
      t267 = epinv * t93
      t270 = (t52 * (t27 * (t339 * t4 + t138) + t255) - t264) * t157
      t271 = t30 * t2 * t52
      t291 = -t99 + t188
      t296 = t291 * t24
      t299 = t28 * t31
      t300 = t26 * t8
      t302 = ddreal(9) * t131
      t304 = t295 * propW16
      t308 = ddreal(1) + t304
      t309 = t40 * propW16
      t311 = t309 * s16
      t312 = t16 + t11 + t311
      t261 = epinv2 * (t20 * (t309 * t33 - t10 + t11) + t22 * (t261 * pr
     &opW16 - ddreal(1)))
      t313 = t1 * t312
      t314 = (t410 + t380) * epinv
      t312 = t94 * ((t20 * t312 + t22 * t308) * epinv - t311 * t254)
      t315 = propW16 * s16
      t316 = t315 * t254
      t317 = t316 * t52 * t244
      t240 = -t240 + t477
      t322 = t25 * t142
      t324 = t322 - t130
      t325 = t25 * t141
      t326 = t325 - t127
      t154 = t25 * t154
      t328 = t154 - t111
      t329 = t153 * t25
      t330 = (t32 - t213) * t157
      t276 = t52 * (t31 * (t328 * s34 + t122 * t445 - t139 * (t149 + t33
     &9) + t225 * t324 + t276 * t326 + t47 * (t329 - t110) - t28) + s34 
     &* (t157 * t240 + t117) * t24 - t194 * t27) + t330 * t24 + t26 * (t
     &147 * t52 - t13 - t14 - t9) + (t101 + t103 + t15) * t24 * epinv
      t308 = t1 * t308
      t331 = t131 * (t308 + t2 + t7 + t8 - t231 - t482)
      t161 = t16 * t276 + t18 * t331 - t40 * (-t192 + ((-t234 - t66 + t4
     &75 - t450) * s16 + t241 + t417) * t52 * t254 * propW16 + t131 * (t
     &147 - t191)) - ddreal(12) * t317 + t157 * (t495 * (t158 - t292 + t
     &419) + t293) + t52 * (t20 * (-t7 * t12 + t12 * (t482 + t231)) + t2
     &4 * (t67 * (epinv * t47 + s34) - t286) + t314 * t27 + t312) + t52 
     &* (t157 * (t261 + t281 + t282) + t20 * (-t8 * t12 + t313) + t24 * 
     &(-t25 * t89 + t47 * (t73 * t86 + t96 + t98)) + t27 * ((t284 + t100
     &) * (s25 + s26 - s34 + s56) + t283) + t31 * (t124 * t445 - t373 * 
     &t81 - t59 * (t149 + t339 + s26 + s56) - t76 * (epinv * t373 + s25 
     &+ s26 + s56) + t339 * t48)) + t204 * ((t71 + t163) * s34 + t161 * 
     &t335 + t165) + t271 + t270 + t204 * (-t92 + t492 - t267 + t175 - t
     &176)
      t276 = t16 + epinv
      t281 = t10 * s56
      t282 = t397 + t21 + t281
      t284 = -t16 * (s16 - s34 + s13 - s45 - s56) + s25 - t69
      t286 = t10 * t359
      t292 = t41 * s26
      t293 = t315 * t10
      t331 = ddreal(1) + t293
      t315 = ddreal(12) * t315
      t333 = t16 + t315 + t11
      t327 = t16 * (s12 + s26 + s16) - t327
      t338 = s12 + s26 + s16 - s34 + s56
      t339 = t10 * t338
      t338 = t16 * t338
      t269 = t398 + t269 + t42 - t386
      t341 = -s13 + s25 + s26 + s45 + s56
      t343 = t385 * t41 + t128 + t286 + t78
      t344 = t339 + s25
      t346 = t25 * t155
      t196 = (-t196 - t197) * t25
      t197 = t15 * t276
      t255 = t52 * (t27 * (t377 * t466 - t138) - t255) + t264
      t264 = t204 * (t216 + t114 + t92 - t492 + t267 - t175)
      t140 = t157 * t255 + t268 + t157 * (t27 * (t205 * t150 - t64 - t65
     &) + t495 * (t193 * t263 - t158 - t419) - t261 * t52) + t20 * (t12 
     &* (-t52 * (t482 + t231 + t51) + t13 + t9) - t313 * t52) + t52 * (t
     &24 * (-t100 * t359 - t101 * (t68 - t339) - t103 * (t266 - t338) + 
     &t157 * (-t257 - t252) + (t346 - t115) * (-t16 * t287 + s25) + (t22
     &8 - t106) * (-t338 + s25)) + t27 * (-t119 * t327 - t132 * t343 - t
     &134 * t269 + t157 * (-t273 - t272) + t25 * (t143 * t343 + t144 * t
     &269 + t152 * t344 + t156 * t327 + t140)) + t31 * ((epinv * t323 + 
     &t421) * s34 + t47 * (epinv * t464 + t434))) + t204 * ((-t164 - t39
     &) * epinv + t139 * (t70 * t341 + s25) + t196) + t294 * (-t108 * t3
     &44 - t283 - t314) - t197 * t24 - t312 * t52 + t264
      t150 = -t172 + t82
      t164 = epinv * t459
      t228 = epinv * t409
      t252 = t173 * t382
      t255 = t97 * (t96 + t98)
      t83 = t27 * (-t205 * t83 * t374 * t157 + t54) + t53 * t6 + t271 + 
     &t52 * (epinv * (t27 * t84 + t465 * (t20 * t333 + t22 * t331)) + t2
     &0 * (t12 * t3 + t224 * t29 + t333 * t450) + t24 * (t67 * (t16 * (e
     &pinv * s26 - s25) + t266 + t407 - t78 - t286 - t292) + t73 * (s25 
     &* t86 - t10 * (s12 + s16 - s34) - t16 * t352) + t85 * (t16 * (-t39
     &3 * t86 + s26 - s34 + s56) + s25) + t181 + t252) + t31 * (-t124 * 
     &t282 + t25 * t484 + t335 * (-t164 - t438) - t59 * (epinv * t282 + 
     &s25 + t74) - t76 * (epinv * t284 + s25 + t74) - t81 * t284 + (t162
     & + t71) * (t16 * (s26 - s13 + s45 + s56) + t21) + t228) + t27 * t1
     &50 * (t342 * t385 + t359 * t40 + s25 + t78) + t27 * (t374 * t5 + t
     &379 * t4) * t157) + t255 * t495 + t204 * (t25 * t489 + t166 - t236
     & - t289 - t423)
      t162 = (-t234 + t475) * s16 + t241 + t417
      t205 = t7 + t8
      t83 = t10 * (t52 * (t195 * t27 + t31 * (t120 * t324 - t122 * t341 
     &+ t285 * t326)) + (t240 * t52 * s34 - t213 + t32) * t157 * t24) - 
     &t16 * t140 + t18 * (t52 * (-t26 * t205 - t299) + t131 * (t331 * t4
     &50 + t2 + t224 - t3)) - ddreal(24) * t317 - t40 * (t296 - t192 + t
     &495 * (-t117 - t118 + t182) * s34 + t131 * (t482 + t231 + t51 - t3
     &08 - t191) - t294 * t186 * (t338 + s25)) + t83 - ddreal(12) * t52 
     &* t254 * propW16 * t162 + ddreal(9) * t131 * t205
      t117 = s25 - s26 + s56
      t140 = t57 + t402 + t42 + s56
      t186 = t86 * s12
      t195 = t16 * (t186 + s25 + s26)
      t234 = s16 + s56
      t61 = t61 + s25 - s13 + s45 + s46
      t240 = t285 + t78 - t77
      t241 = epinv * t50
      t77 = t204 * ((t409 - t39) * epinv + t459 * (-t57 * t25 - s16 - s5
     &6 - t21 + t266 - t292) - t76 * (epinv * (t385 - t78 + t77) + t195)
     & - t81 * (t385 - t78 + t259) + t224 + t28 + t332 + t337 + t92)
      t50 = t50 * t71
      t23 = t52 * (t20 * ((t4 + t5) * t12 - t180) + t24 * (-t96 * t97 - 
     &t178 + t179) + t31 * ((t93 - t94) * epinv - t492 - t50)) + t157 * 
     &(t27 * (-t169 + t168 - t64 - t65) - t170) + t24 * (-t171 * t52 - t
     &54 * t86 - t449) + t53 * (-t3 - t51) + t52 * (t20 * (t12 * t231 - 
     &t12 * t250 - t29 * t423) + t27 * (-t103 * (t265 + t21 + t87) - t13
     & * (t57 + t33) - t256 * t98 - t382 * t90 - t73 * (epinv * t280 + t
     &195) - t88 * t91 + (t172 - t82) * (t57 + s26 - s16 + s34 - s56 + t
     &21) - t494) + t31 * (-t117 * t63 - t140 * t72 - t17 * t489 + t25 *
     & (t140 * t34 - t468 - t95) + t438 * t58 + t46 * (s26 * t276 + t149
     & - t353) - t48 * (t241 - s25) - t484 * t91 + t59 * (t16 * (epinv *
     & t61 - s26) - t21)) - t288 * t24) - t23 * t243 + t77
      t29 = t25 * t156
      t77 = t29 - t119
      t86 = t25 * t357
      t87 = t323 * s34
      t140 = t25 * t143
      t143 = (t140 - t132) * t133
      t156 = t434 * t31
      t152 = t25 * t152
      t112 = t52 * (t27 * (-t190 - t189) + t31 * ((t230 - t418) * s34 - 
     &t184 * t321 + t97 * (t112 - t345))) + t131 * (-s25 * t9 - t2 - t23
     &6 + t7)
      t15 = -t10 * t52 * ((t152 - t108) * t35 * t27 + t31 * t285 * (t329
     & - t110)) - t16 * (t20 * (t12 * t9 + t52 * (-t7 * t12 - t217)) + t
     &24 * (t52 * (-t210 - t212 + t211 - t177) - t330 - t201 * t157 * t5
     &2 - t202 * t157 * t52) + t27 * (t52 * (t135 * (t144 * t25 - t134) 
     &- t206 + t143) - t106 - t115 - t15) + t27 * (t25 * (t151 + t155) +
     & t77 * t52 * t225) - t30 * t236 * t52 + t204 * (t117 * (t25 * (t14
     &1 + t142) - t127 - t130) + t122 * t240 + t124 * t61 + t139 * (epin
     &v * t240 - s25 + s26 - s56) + t225 * t328 - t25 * t416 + t464 * (t
     &68 - t37) + (t86 - t480) * (t57 + t234 + t21 + t292) + (t25 * t383
     & - t433) * (s16 + s34 + s56 + t42 + t36) + t187 - t87) + t156) + t
     &18 * (t26 * t6 * t52 + t131 * (t423 + t243 - t250 + t231 + t1 + t4
     & + t5)) - t40 * t112 - t23 - t302 * t6
      t23 = ddreal(1) - t293
      t30 = t16 - t315 + t11
      t35 = ddreal(1) + t70
      t37 = -s13 + s25 + s45 + s46 + s56
      t42 = t16 * (s25 - s12 - s16 + s34) + s56 - t42
      t61 = -t16 * t425 + s25
      t68 = -s12 - s13 + s25 + s45 + s46 + s56
      t106 = s12 - s13 + s25 + s45 + s46 + s56
      t112 = t16 * t47
      t134 = t402 + t78 - s26
      t44 = t74 + t44 + t281
      t135 = s25 * t40 + t10 * (s34 + s56) - t369 - t57
      t144 = t47 + t295
      t151 = t10 * t385 + t16 * t234 + s25
      t11 = t496 * (t27 * (-t52 * (t208 + t250 + t207) - t197) + t31 * (
     &t490 * (-t246 - t249) + t52 * (-t91 * (t484 + t489) - t2 + t224 + 
     &t228 - t28)) + t222 * (t20 * (t16 + t11 - t311) + t22 * (ddreal(1)
     & - t304)) - epinv2 * (t20 * (t309 * t47 + t11 + t16) + t22 * (t112
     & * propW16 + ddreal(1))) * t52 * t490 + t495 * (-t101 * (epinv * t
     &262 - t128 - t69) + t211 - t231 - t288 - t482 - t497 - t51))
      t9 = t298 * (t24 * (-t13 * t33 - t262 * t9 - t218) + t31 * (-t110 
     &* t44 - t127 * t134 - t130 * t42 - t139 * (epinv * t123 - s26 + s5
     &6 + t21) + t490 * (t413 + t414 + t415)) - t103 * t27 * (-t265 + t3
     &6))
      t9 = t9 + t11 + t496 * (t31 * (t490 * (t52 * (t25 * (-t335 * t468 
     &- t290 + t306) - t337 * t340 - t412) + t491) - t52 * (t215 + t214)
     & + t52 * (t154 - t111) * (-t369 - t113 + t386 + s25)) + t294 * (t1
     &21 * (t29 - t119) + (t140 - t132) * t133 - t107 - t209 - t494 + t1
     &09 * t152) + t495 * (-t100 * (s12 - s16 - s26 + s34 - s56) + t301 
     &* (t106 * t56 + t136 * t68 + t37 * t55))) + t496 * (t24 * (t52 * (
     &-t157 * (t106 * t147 + t146 * t37 + t32 * t68) - t283 * epinv) - t
     &449) + t27 * (t52 * (t25 * (t200 + t198 + t199) - t206 * t35) - t1
     &4) + t204 * (t25 * (t134 * t141 + t142 * t42 + t153 * t44 + t274 *
     & t490) + t490 * (t332 * t335 + t275 + t277 - t307) + t279 * t490))
      t11 = (t157)**2
      t13 = t496 * (t24 * (t35 * t54 + t52 * (t137 * t157 + t151 * t96 +
     & t67 * (epinv * t151 - s26 - s34 + s56 - t128 - t69) + t85 * (-t16
     & * (-t186 + t265) + s25) - t179 + t4 + t5)) + t52 * (-t1 * (t20 * 
     &t30 + t22 * t23) + t94 * ((-t20 * t30 - t22 * t23) * epinv - ddrea
     &l(12) * t316)) + t204 * (t135 * t72 + t144 * t63 + t46 * (epinv * 
     &t144 - t112) + t48 * (t241 - t128) + t59 * (t16 * (-epinv * t62 + 
     &s26) + t281 + t371) + t176 + t492 + t50 - t92))
      t22 = t496 * (t27 * (t52 * (t185 - t105) * t58 - t99 + t188 + t52 
     &* (t182 - t118) * s34) + t204 * (t183 - t125) * t117)
      t23 = t298 * (t31 * (t219 * t47 + t233) + t254 * t162 * propW16)
      t29 = t464 * t47 + t87
      t30 = t298 * (t174 * t27 + t24 * (t177 + t3))
      t9 = -t10 * t22 + t16 * t30 + t258 * (t13 + t298 * (t27 * (t150 * 
     &t45 + t157 * (-t148 * (t25 * t37 + t186) + (-t168 + t64 + t65) * t
     &68) + t61 * t98 + t73 * (epinv * t61 + t21 + t69 + t74) + t90 * (e
     &pinv * t10 - t41) + t252 + (s12 + s13 + s16 + s26 - s34 - s45 - s4
     &6) * (-t159 + t138) * t11) + t31 * (t76 * (-t10 * t238 + t21 - t74
     & + t80) - t267 + t79 - t160 * t135 + t39 * t17) + t158 * t24 * t15
     &7)) - ddreal(4) / ddreal(3) * t9 - ddreal(8) / ddreal(3) * t496 * 
     &(t31 * (t490 * (-t25 * t411 + t248 + (-t227 + t474) * t52 * s34) +
     & t52 * (t25 * (t334 * t383 - t318 - t416) - t334 * t433 - t335 * t
     &459 + t187 - t236 - t423 + t483) + t52 * (t86 - t480) * (t112 + s3
     &4)) - t294 * (s34 * t408 + t6)) - ddreal(16) * t316 * t298 * t244 
     &+ ddreal(16) / ddreal(3) * t496 * t204 * t29 - ddreal(8) * t23
      t13 = t319 + t69 + s25
      t22 = epinv * t97
      t23 = t319 + s25
      t30 = -t17
      t35 = t6 - t51
      t13 = t24 * (-t73 * (t74 + t69 + t21 - t22) - t90 * t91 + (t88 + t
     &89) * t30 - t179 + t255 + t3 + t4 + t5) + t31 * (-t38 * t72 + t2 -
     & t224 - t228 + t236 - t28 + t423 - t450) + t300 + t24 * (s25 * t54
     & + t45 * (t172 - t82) - t67 * (-t22 + s26 + s34 - s56 + t69 + t128
     &) + t177 + t178) + t31 * (t124 * t23 + t13 * t81 + t25 * (t34 * t3
     &8 - t484 - t489) + t335 * (t164 + t438) - t59 * (-epinv * t23 + s2
     &5 + t74) - t76 * (-epinv * t13 + s25 + t74) + t95 * t30 + t167 - t
     &66) + t26 * t7
      t11 = struc9MP * (t10 * t31 * (t120 * (t322 - t130) - t122 * t49 +
     & t285 * (t325 - t127)) + t16 * (t24 * (t212 - t107 - t209 - t211 +
     & t191) + t24 * (t101 * t102 + t103 * t104 + t109 * (t152 - t108) +
     & t116 * (t346 - t115) + t121 * t77 + t25 * (t198 + t200) + t143 + 
     &t203 * s25) + t31 * (epinv * t29 + s34 * t421 - t139 * (t70 * t49 
     &+ s25) - t114 - t196 - t216) + t206 * t27 + t26 * t35 + t156 * t47
     &) - t18 * t253 * t205 + t40 * (t24 * (s25 * t291 - t189 - t190) - 
     &t253 * t35 + t192 * s25) + t13) + struc60MP * t27 * (-t33 * (-t16 
     &* (t146 + t147 + t32 - t194 - t213 - t193) - t148 * t25 - t168 + t
     &64 + t65) + t138 - t159) * t11
      t3 = t496 * (t52 * t11 + (t16 * t129 + t18 * t131 * (t1 + t3 + t4 
     &+ t5 - t7 - t8) + t40 * (t27 * (t189 * t52 + t190 * t52 - t188 + t
     &99) + t192 * t52 * t126 + (-t6 + t2 + t191) * propZ25 * t19) + t52
     & * (t20 * (t12 * t174 - (t3 + t4 + t5) * t12 + t180) + t24 * (t97 
     &* (t96 + t98) + t177 + t178 - t179) + t27 * (t17 * t173 + t45 * t8
     &2 + t181) + t31 * ((-t93 + t94) * epinv - t59 * (t16 * (epinv * t6
     &2 + s26) + s25) + t79 + t175 - t176)) + t157 * (t27 * (t169 - t168
     & + t64 + t65) + t170) + t24 * t54 + t51 * t53 + t52 * (t24 * (t73 
     &* ((-t18 + epinv) * s25 - t69 - t74 * t25) + t171) - t172 * t27 * 
     &t45) + t52 * ((t75 + t163) * s34 + t76 * (-t74 + t80 - s25) - t92 
     &+ t165 + t166 + t167) * t31) * struc10MP + t83 * struc5MP + t15 * 
     &struc6MP)
      t4 = (-t10 * t298 * t451 * t254 * t244 + t16 * t242 + t245 * (t145
     & + t43 + t298 * t297 + t496 * (t31 * (t490 * (t52 * (t247 * (t412 
     &- t413 - t414 - t415) + t332 * t336 + t351 * (-t16 * (s13 * (-t149
     & - s26 - s56) + t278) - s25 * (t149 + t16 * (s45 * t25 + s46 * t25
     &) + t78 + t352 - t353))) - t491 * t247) - t52 * (t160 * t310 + t49
     &2 * t223) + t52 * (t223 * t93 - t229 * t409) * epinv) + t52 * (t15
     &7 * (t305 * (t64 + t65) + t376 * t5 + t4 * (s13 * t379 + t220 * t3
     &78 + t57 * t378) - t138 * t303 * t157) - t172 * t237 + t493 * t320
     & + t494 * t221) * t27 + t495 * ((t226 * t84 + t235 * t410) * epinv
     & + t157 * (t137 * t420 + t260 * t419)))) + t258 * t60 + t298 * t8 
     &* t27 * t320) * struc1MP + t9 * struc7MP
      result = t258 * t3 + ddreal(4) / ddreal(3) * t496 * ((-t16 * t232 
     &+ t18 * (t52 * (-t299 - t300) + t131 * (t450 + t243 - t231 + t7 + 
     &t2)) - t40 * (t296 - t192 + t131 * (-s25 * t14 - t1 + t147 + t482 
     &- t483)) + t239 + t271 + t270 + t204 * (s34 * (t75 + t25 * t490 * 
     &(-t468 + t251)) - t92 - t267) + t302 * t8) * struc2MP + (struc3MP 
     &+ struc8MP) * t161 + t204 * (-t16 * (epinv * t139 + t122) - (t59 +
     & t76 - t48) * epinv - t124 + t71 - t81) * struc61MP) + t4

           ampNonresonantLightFullImC4MP = result
       end function ampNonresonantLightFullImC4MP

       function ampNonresonantLightFullImC4PM()
           implicit none
           type(dd_complex) :: ampNonresonantLightFullImC4PM
           type(dd_complex) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           type(dd_complex) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           type(dd_complex) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           type(dd_complex) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           type(dd_complex) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           type(dd_complex) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           type(dd_complex) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           type(dd_complex) ::  t175,t176,t177,t178,t179,t18,t180,t181,t182,t183,t184,t185
           type(dd_complex) ::  t186,t187,t188,t189,t19,t190,t191,t192,t193,t194,t195,t196
           type(dd_complex) ::  t197,t198,t199,t2,t20,t200,t201,t202,t203,t204,t205,t206
           type(dd_complex) ::  t207,t208,t209,t21,t210,t211,t212,t213,t214,t215,t216,t217
           type(dd_complex) ::  t218,t219,t22,t220,t221,t222,t223,t224,t225,t226,t227,t228
           type(dd_complex) ::  t229,t23,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239
           type(dd_complex) ::  t24,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249,t25
           type(dd_complex) ::  t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t26,t260
           type(dd_complex) ::  t261,t262,t263,t264,t265,t266,t267,t268,t269,t27,t270,t271
           type(dd_complex) ::  t272,t273,t274,t275,t276,t277,t278,t279,t28,t280,t281,t282
           type(dd_complex) ::  t283,t284,t285,t286,t287,t288,t289,t29,t290,t291,t292,t293
           type(dd_complex) ::  t294,t295,t296,t297,t298,t299,t3,t30,t300,t301,t302,t303
           type(dd_complex) ::  t304,t305,t306,t307,t308,t309,t31,t310,t311,t312,t313,t314
           type(dd_complex) ::  t315,t316,t317,t318,t319,t32,t320,t321,t322,t323,t324,t325
           type(dd_complex) ::  t326,t327,t328,t329,t33,t330,t331,t332,t333,t334,t335,t336
           type(dd_complex) ::  t337,t338,t339,t34,t340,t341,t342,t343,t344,t345,t346,t347
           type(dd_complex) ::  t348,t349,t35,t350,t351,t352,t353,t354,t355,t356,t357,t358
           type(dd_complex) ::  t359,t36,t360,t361,t362,t363,t364,t365,t366,t367,t368,t369
           type(dd_complex) ::  t37,t370,t371,t372,t373,t374,t375,t376,t377,t378,t379,t38
           type(dd_complex) ::  t380,t381,t382,t383,t384,t385,t386,t387,t388,t389,t39,t390
           type(dd_complex) ::  t391,t392,t393,t394,t395,t396,t397,t398,t399,t4,t40,t400
           type(dd_complex) ::  t401,t402,t403,t404,t405,t406,t407,t408,t409,t41,t410,t411
           type(dd_complex) ::  t412,t413,t414,t415,t416,t417,t418,t419,t42,t420,t421,t422
           type(dd_complex) ::  t423,t424,t425,t426,t427,t428,t429,t43,t430,t431,t432,t433
           type(dd_complex) ::  t434,t435,t436,t437,t438,t439,t44,t440,t441,t442,t443,t444
           type(dd_complex) ::  t445,t446,t447,t448,t449,t45,t450,t451,t452,t453,t454,t455
           type(dd_complex) ::  t456,t457,t458,t459,t46,t460,t461,t462,t463,t464,t465,t466
           type(dd_complex) ::  t467,t468,t469,t47,t470,t471,t472,t473,t474,t475,t476,t477
           type(dd_complex) ::  t478,t479,t48,t480,t481,t482,t483,t484,t485,t486,t487,t49
           type(dd_complex) ::  t5,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t6
           type(dd_complex) ::  t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t7,t70
           type(dd_complex) ::  t71,t72,t73,t74,t75,t76,t77,t78,t79,t8,t80,t81
           type(dd_complex) ::  t82,t83,t84,t85,t86,t87,t88,t89,t9,t90,t91,t92
           type(dd_complex) ::  t93,t94,t95,t96,t97,t98,t99

           type(dd_complex) :: result

      t1 = intHLs160000x0111D0eps0()
      t2 = intHLs16s25s26s34s56x1112D4eps0()
      t3 = intHs160000x0112D2eps0()
      t4 = intHs160s26s34s56x1021D2eps1()
      t5 = intHs160s26s34s56x1022D4eps0()
      t6 = intHs160s26s34s56x1031D4eps0()
      t7 = intHs16s25s26s34s56x1121D4eps0()
      t8 = intHs16s25s26s34s56x1211D4eps0()
      t9 = ddreal(1) - epinv
      t10 = ddreal(4)
      t11 = propZ25 * s25
      t12 = -t11 + t10
      t13 = ddreal(3)
      t14 = (gw)**2
      t15 = (gb)**2
      t16 = t13 * s25
      t17 = t16 * t14 * propZ25
      t18 = t15 * t12 + t17
      t19 = intHs16s25s26s34s56x1311D4eps1()
      t20 = ddreal(2)
      t21 = -t20 + epinv
      t22 = -t12
      t23 = t15 * t22
      t24 = -t17 + t23
      t25 = intHLs16s25s26s34s56x1111D2eps0()
      t26 = t11 + t20
      t27 = t15 * t26
      t28 = t17 + t27
      t29 = intHs16s25s26s34s56x1310D4eps0()
      t30 = s16 + s26 - s34 + s56
      t31 = s12 + s25
      t32 = t20 * t31
      t33 = t10 * s26
      t34 = t33 + t32 + s16
      t35 = intHLs16s25s26s34s56x1211D2eps1()
      t36 = ddreal(6)
      t37 = ddreal(5)
      t38 = t13 * s26
      t39 = t37 * s25
      t40 = t36 * s12
      t41 = t40 + t39 + t38 + s34 - s56
      t42 = intHLs16s25s26s34s56x1121D2eps1()
      t43 = s25 + s26 + s56
      t44 = intHLs16s25s26s34s56x1112D2eps1()
      t45 = -s13 + s26 + s45 + s46 + s56
      t46 = t20 * t45 + s25
      t47 = intHs16s25s26s34s56x1220D4eps1()
      t48 = intHs16s25s26s34s56x1130D4eps1()
      t49 = t20 * s12
      t50 = s25 + t49
      t51 = intHs16s25s26s34s56x1112D4eps0()
      t52 = ddreal(1) / s25
      t53 = t15 * (-propZ25 * t20 + ddreal(8) * t52) + t36 * t14 * propZ
     &25
      t54 = intHs16s25s26s34s56x1122D4eps1()
      t55 = t21 * s25
      t56 = t10 * s12
      t57 = t20 * epinv
      t58 = intHLs16s25s26s34s56x1112D2eps0()
      t59 = intHLs16s25s26s34s56x1122D4eps1()
      t60 = s13 - s45 - s46
      t61 = t20 * s56
      t62 = -t60 + t61
      t63 = intHs16s25s26s34s56x1211D2eps1()
      t64 = intHs16s25s26s34s56x1321D6eps1()
      t65 = t16 + t49
      t66 = intHs16s25s26s34s56x1231D6eps1()
      t67 = t56 + t16
      t68 = intHLs16s25s26s34s56x1132D6eps1()
      t69 = t20 * s26
      t70 = t69 + s25
      t71 = intHLs16s25s26s34s56x1312D6eps1()
      t72 = s25 + s16
      t73 = t20 * (s12 + s26)
      t74 = t72 + t73
      t75 = intHs16s25s26s34s56x1312D6eps1()
      t76 = t16 + t73
      t77 = intHs16s25s26s34s56x1213D6eps1()
      t78 = s12 + s25 + s26
      t79 = t20 * t78 + s34
      t80 = s25 + s26 - s56
      t81 = intHLs16s25s26s34s56x1123D6eps1()
      t82 = t20 * s25
      t83 = t38 + t82 - s56
      t84 = intHLs16s25s26s34s56x1213D6eps1()
      t85 = t38 + t32 - s56
      t86 = intHs16s25s26s34s56x1132D6eps1()
      t87 = t56 + s26 + s34 - s56 + t82
      t88 = intHs16s25s26s34s56x1123D6eps1()
      t89 = t20 * (s12 + s34) + t80
      t90 = intHs160s26s34s56x1020D2eps1()
      t91 = intHs16s25s26s34s56x1110D2eps1()
      t92 = intHs16s25s26s34s56x1220D4eps0()
      t93 = intHs16s25s26s34s56x1130D4eps0()
      t94 = intHs16s25s26s34s56x1120D2eps1()
      t95 = t9 * s25
      t96 = intHLs16s25s26s34s56x1212D4eps0()
      t97 = s25 + s26 - s34
      t98 = t20 * t60
      t99 = t13 * s56
      t100 = t56 + t98 - t99 + t97
      t101 = intHs16s25s26s34s56x1222D6eps0()
      t102 = intHLs16s25s26s34s56x1121D2eps0()
      t103 = intHs16s25s26s34s56x1120D2eps0()
      t104 = intHs16s25s26s34s56x1210D2eps0()
      t105 = intHs16s25s26s34s56x1210D2eps1()
      t106 = intHs160000x0112D2eps1()
      t107 = intHs16s25s26s34s56x1121D2eps1()
      t108 = ddreal(1) + epinv
      t109 = t49 * t108
      t110 = intHLs160000x0112D2eps0()
      t111 = intHLs160000x0111D0eps1()
      t112 = intHLs160000x0112D2eps1()
      t113 = intHLs16s25s26s34s56x1112D4eps1()
      t114 = -t13 + t57
      t115 = intHs16s25s26s34s56x1122D4eps0()
      t116 = intHs16s25s26s34s56x1212D4eps0()
      t117 = intHLs16s25s26s34s56x1222D6eps0()
      t118 = intHs16s25s26s34s56x1212D4eps1()
      t119 = intHLs16s25s26s34s56x1211D2eps0()
      t120 = intHLs16s25s26s34s56x1212D4eps1()
      t121 = intHs16s25s26s34s56x1121D4eps1()
      t122 = intHs16s25s26s34s56x1211D4eps1()
      t123 = intHs16s25s26s34s56x1112D4eps1()
      t124 = intHs160s26s34s56x1022D4eps1()
      t125 = intHs16s25s26s34s56x1411D6eps0()
      t126 = intHs16s25s26s34s56x1131D4eps0()
      t127 = intHs16s25s26s34s56x1221D4eps0()
      t128 = intHs16s25s26s34s56x1311D4eps0()
      t129 = intHs16s25s26s34s56x1141D6eps0()
      t130 = intHs16s25s26s34s56x1321D6eps0()
      t131 = intHs16s25s26s34s56x1231D6eps0()
      t132 = intHs160s26s34s56x1031D4eps1()
      t133 = intHs16s25s26s34s56x1121D2eps0()
      t134 = intHs16s25s26s34s56x1131D4eps1()
      t135 = t55 - t56
      t136 = intHs16s25s26s34s56x1221D4eps1()
      t137 = t55 - t49
      t138 = intHLs16s25s26s34s56x1132D6eps0()
      t139 = intHLs16s25s26s34s56x1312D6eps0()
      t140 = intHs16s25s26s34s56x1312D6eps0()
      t141 = intHs16s25s26s34s56x1113D4eps1()
      t142 = intHs16s25s26s34s56x1114D6eps0()
      t143 = intHs16s25s26s34s56x1213D6eps0()
      t144 = intHLs16s25s26s34s56x1113D4eps0()
      t145 = t98 - t99 - s26
      t146 = intHLs16s25s26s34s56x1122D4eps0()
      t147 = intHLs16s25s26s34s56x1114D6eps0()
      t148 = intHLs16s25s26s34s56x1123D6eps0()
      t149 = intHLs16s25s26s34s56x1213D6eps0()
      t150 = intHs16s25s26s34s56x1132D6eps0()
      t151 = intHs16s25s26s34s56x1123D6eps0()
      t152 = intHs160s26s34s56x1020D2eps0()
      t153 = intHs16s25s26s34s56x1310D4eps1()
      t154 = intHs16s25s26s34s56x1110D2eps0()
      t155 = intHLs16s25s26s34s56x1113D4eps1()
      t156 = ddreal(1) / t30
      t157 = t21 * t91
      t158 = t94 * (t95 + t49) + t157
      t159 = t9 * intHs16s25s26s34s56x1222D6eps1()
      t160 = (-epinv * t44 - t58) * t46
      t161 = t35 * s34
      t162 = t9 * intHLs16s25s26s34s56x1222D6eps1()
      t163 = epinv * t106
      t164 = t107 * (t109 + s25)
      t165 = t54 * (t57 * s26 - s26 - s34 + s56 + t55 - t56)
      t166 = (-t9 * (t52 * t90 + t105) + t103 + t104) * t18
      t167 = (-t159 + t101) * t41
      t168 = t3 + t4 + t5
      t169 = t121 + t122
      t170 = t7 + t8
      t171 = t113 * t114
      t172 = s34 * t119
      t173 = t9 * t124
      t174 = t154 * t156
      t175 = t152 * t156
      t176 = t123 * t114
      t177 = t1 * t26
      t178 = t9 * intHLs16s25s26s34s56x1114D6eps1()
      t179 = -t178 + t147
      t180 = t9 * intHs16s25s26s34s56x1141D6eps1()
      t181 = -t180 + t129
      t182 = t9 * intHs16s25s26s34s56x1114D6eps1()
      t183 = intHs16s25s26s34s56x1411D6eps1() * t9
      t184 = (-t182 + t142) * s34
      t185 = t181 * t50
      t186 = t179 * t28
      t187 = t64 * t65
      t188 = t66 * t67
      t189 = t75 * t76
      t190 = t77 * t79
      t191 = t86 * t87
      t192 = t88 * t89
      t193 = t68 * t70
      t194 = t71 * t74
      t195 = t19 * t21
      t196 = t2 * t26
      t197 = t130 * t65
      t198 = s12 * t133
      t199 = t134 * t135
      t200 = t9 * t153
      t201 = t70 * t138
      t74 = t74 * t139
      t202 = t131 * t67
      t203 = s34 * t141
      t76 = t76 * t140
      t79 = t79 * t143
      t204 = t87 * t150
      t205 = t89 * t151
      t206 = t52 * t18
      t83 = t24 * (t52 * (t136 * t137 + t9 * (-t156 * (t31 * t47 + t48 *
     & t50) - t132) + t198 + t199) - t200 * t156) + t24 * (t156 * (t52 *
     & (t31 * t92 + t50 * t93) + t29) + t126 + t127 + t128 + t195) + t52
     & * (t15 * (t22 * t6 + t196) - t25 * t28 + t197 * t18) + t52 * (t18
     & * (-t187 - t188 - t189 - t190 - t191 - t192) + t28 * (-t81 * t83 
     &- t84 * t85 - t193 - t194)) * t9 + t206 * (t79 + t204 + t205 + t20
     &2 + t203 + t76) + (-t144 * t145 + t146 * t62 + t148 * t83 + t149 *
     & t85 + t155 * (-epinv * t145 + s25 + s26 - s56) + t201 + t74) * t5
     &2 * t28
      t85 = t14 * propZ25
      t14 = t13 * t85 * (t1 - t3 - t4 - t5 + t7 + t8) + t20 * t83 + t36 
     &* (t18 * (t184 * t52 + t185 * t52 + t125 - t183) + t186 * t52 * t8
     &0 + (t2 - t6) * propZ25 * t14) + t52 * (t15 * (t12 * t170 + t168 *
     & t22 + t177) + t18 * (t169 * t21 + t175 + t176) + t24 * (t118 * ((
     &-t13 + epinv) * s25 - t56 - t69 * t9) + t70 * (t115 + t116) - t173
     & + t174) + t28 * ((t111 - t112) * epinv + t117 * t34 - t120 * (epi
     &nv * t100 - s25 - t69) + t171 - t172)) + t156 * (t52 * t158 * t24 
     &+ t166) + t24 * (t52 * (t165 + t164 + t163) + t63) + t51 * t53 + t
     &167 * t52 * t18 + t52 * (t43 * (-epinv * t42 - t102) + t59 * (t20 
     &* (epinv * t62 + s26) + s25) - t96 * t100 - t110 + t160 - t161 * e
     &pinv - t162 * t34) * t28
      t62 = intHL0s25s26s34s56x1130D4eps0()
      t80 = -s13 + s25 + s26 + s56
      t83 = s45 + s46
      t100 = t20 * t83
      t145 = t80 + t100
      t207 = s13 - s45 - s56
      t208 = -t20 * t207 + s25 + s46
      t209 = intHLs160000x0211D2eps0()
      t210 = -s16 + s25 + s34
      t211 = t20 * (s12 + s13 - s56)
      t212 = t13 * s46
      t213 = t10 * s45
      t214 = t213 + t212 + t210 - t211
      t215 = t10 * t83
      t211 = t215 + t210 - t211
      t216 = intHs160000x0121D2eps0()
      t217 = s25 - s26 - s16 + s34
      t218 = s13 - s56
      t219 = t13 * t218
      t220 = t215 - t219 + t217 - t49
      t221 = intHLs16s25s26s34s56x1211D4eps0()
      t222 = s25 - s26
      t223 = s16 - s34
      t224 = t20 * t223
      t225 = t37 * s46
      t226 = t36 * s45
      t227 = -t56 + t225 - t219 + t222 - t224 + t226
      t228 = intHs16s25s26s34s56x1111D2eps0()
      t229 = s26 + s16 - s34 + s13 - s56
      t230 = t20 * (s12 - s45 - s46)
      t231 = t229 + t230
      t232 = intHL0s25s260s56x1012D2eps0()
      t233 = -t98 + t43
      t234 = intHL0s25s260s56x1013D4eps0()
      t235 = intHL0s25s260s56x1021D2eps0()
      t236 = intHs160000x0111D0eps0()
      t237 = intHL0s25s26s34s56x1210D2eps1()
      t238 = intHL0s25s26s34s56x1310D4eps1()
      t239 = t85 * s25
      t240 = t239 + t15
      t241 = s25 + s16 - s34
      t242 = t241 + t73
      t243 = intHs160s26s34s56x1012D2eps1()
      t244 = s26 + s16 - s34 - s56
      t245 = t20 * (s12 + s13 - s45 - s46)
      t246 = t244 + t245
      t247 = s25 - s12
      t248 = t20 * t247
      t249 = t36 * t83
      t250 = t249 - t219 + s26 - s16 + s34 + t248
      t251 = intHLs16s25s26s34s56x1121D4eps0()
      t252 = s46 + s56
      t253 = t10 * s13
      t226 = t252 * t37 + s26 + t16 + t226 - t253
      t254 = s16 + s26 - s34
      t255 = s56 * t37 + t248 + t249 - t253 - t254
      t225 = -t207 * t36 + t16 + t225
      t256 = s13 + s16
      t257 = t256 + t49
      t258 = -s16 + s25
      t259 = t13 * s13
      t260 = t20 * s34
      t261 = -t260 - t258
      t262 = t13 * s34
      t263 = -t262 - t258
      t264 = t20 * t261 * t83 - (s13 + t49 + s16) * s16 + s25 * t257 + s
     &26 * t258 - (t20 * (s25 - s12 - s16) - s26 - t259 + s34) * s34 + s
     &56 * t263
      t265 = s13 - s45
      t266 = t20 * t265
      t267 = (t257 + s26 - s34) * s25
      t268 = (s46)**2
      t269 = s25 * s45
      t270 = t20 * (-t269 + t268) + (-t266 + s26) * s46 - (-t212 + s25) 
     &* s56 + t267
      t271 = -s13 + s26 + s45
      t272 = s12 + s16
      t273 = t20 * t272
      t274 = s25 - s13
      t275 = t10 * s25
      t276 = t274 * s34
      t277 = s13 * s16
      t278 = t20 * (t277 + t276) - (s26 + t259 + t273 - s25) * s25 + (t2
     &75 - t224) * s45 - (-t275 + t224 + s26) * s46 + (t16 - t224 + s46)
     & * s56
      t279 = s26 - s34
      t280 = t20 * (s16 - s13 + s45) + t279
      t281 = s12 + s16 + s26 - s34
      t282 = t20 * t281
      t283 = (-s13 + s25 + s26) * s34
      t284 = t272 * s13
      t285 = (s12 - s13 + s16 + s25 + s26) * s26
      t286 = (t212 + t282) * s56
      t287 = t20 * (s25 * t272 - s45 * t210 + t268 - t283 - t284 + t285)
      t288 = s16 - s13
      t289 = t37 * s13
      t290 = t289 - t82 - s16
      t291 = s26 + s16
      t292 = s26 + s45 + s46
      t293 = t13 * t288
      t294 = t10 * t292
      t295 = t37 * s34
      t296 = (s56)**2
      t297 = t13 * t296
      t298 = t82 * t288
      t299 = t10 * t284 + (t290 - s26) * s26 + (-t289 + t16 + s26) * s34
     & + (t295 - t294 - t293 - t82) * s56 + t83 * (s34 * t36 - t10 * t29
     &1) - t297 - t298
      t300 = t13 * s16
      t301 = ddreal(7)
      t302 = -s13 * t301 + s26 + t275
      t303 = t301 * s34
      t304 = (s25)**2
      t305 = s12 * s13
      t222 = t20 * (t305 + t304) - t297 - (t300 + s13) * s25 - (-t289 + 
     &t72 + s26) * s26 + t302 * s34 + (-s16 * t37 + s25 + t259 - t294 + 
     &t303) * s56 + t83 * (-s16 * t36 + t10 * t222 + ddreal(8) * s34) + 
     &t289 * s16
      t306 = s12 + s13 + s16
      t307 = s25 + s26
      t308 = t20 * s16
      t309 = t306 * s25
      t310 = (s12 + s16 + s25 + s26) * s26
      t311 = t70 * s34
      t312 = t82 + s26
      t313 = t20 * (-s34 * t307 - s45 * t312 + (-s25 + s12 + s13 - s45 +
     & s46 - s56) * s56 + t268 + t310) + (t56 + t259 + t308) * s25 - (t2
     &66 + s25) * s46
      t314 = intHLs16s25s26s34s56x1321D6eps0()
      t315 = s25 + s16 - s13
      t316 = t10 * s34
      t248 = -t20 * (t83 * (t260 + t307) + t305) + (s16 + t49) * s25 + (
     &t315 + t49 + s26) * s26 + (t259 - t248 + s16 - s34) * s34 + (-t316
     & - s25 + s16 + s13 + t230 - s56) * s56 - t277
      t317 = (-t253 + t16 + t69) * s34
      t318 = t288 * s25
      t292 = -t20 * (t83 * (-t262 + t291) - t284) + (t259 - t72 - s26) *
     & s26 + (-t20 * t292 - s56 - t315 + t316) * s56 + t317 - t318
      t315 = intHLs16s25s26s34s56x1131D4eps0()
      t319 = s12 + s13
      t320 = t20 * t319
      t321 = s12 - s45
      t322 = t20 * t321
      t323 = (t320 + s16) * s25
      t324 = (s25 + s16 + s13 + t49 + s26) * s26
      t325 = t20 * t307
      t326 = intHLs16s25s26s34s56x1221D4eps0()
      t327 = -t38 - t82 + s13
      t328 = s25 + s26 + s34
      t329 = -s25 + s16 - s34 + s13 - s46 + t322
      t330 = t328 * s46
      t331 = t20 * (s45 * t328 - t305)
      t332 = t20 * (s26 - s13)
      t333 = -t332 - s25
      t334 = s16 * s25
      t335 = t20 * ((t20 * t291 + s12 - s13 + s25 - s34 + s45 + s46 + s5
     &6) * s56 + t83 * (t308 + t279) + t285 - t305 + t334)
      t336 = t253 * s16
      t337 = t10 * t319
      t338 = t20 * (s26 + s16 - s13 + s45 + s46)
      t339 = (s34)**2
      t340 = t83 * (-t260 + t291)
      t341 = t20 * (-t277 + t340 - t339 + t296) + (t256 + s26) * s25 + (
     &t337 + t38 + t308) * s34 + (-t262 + t338 + s25) * s56
      t342 = t20 * s13
      t343 = s25 - t342
      t344 = s25 - s34
      t345 = t83 * t254
      t291 = t20 * (-t277 + t345 + t296) + s25 * t256 + s26 * t343 + (t2
     &91 + t320 - s34) * s34 + (t344 + t338) * s56
      t338 = intHLs16s25s26s34s56x1311D4eps0()
      t346 = s45 * t20 - s13 + s25 + s26 + s46 + s56
      t347 = intHLs160000x0111D2eps0()
      t348 = intHLs16s25s26s34s56x1411D6eps0()
      t349 = intHs160s26s34s56x1011D2eps0()
      t350 = intHL0s25s260s56x1022D4eps0()
      t351 = s26 + s56
      t352 = t351 + t82
      t353 = s45 + s46 + s56
      t354 = t20 * t353
      t355 = s13 * t352
      t356 = t351 * t353
      t357 = (-t354 - s26 - s25) * s25 + t355 - t356
      t358 = t241 + t69
      t359 = t319 * s25
      t360 = propW16 * t242
      t361 = t360 * t43 * t30
      t362 = t83 * t210
      t354 = t354 + s25
      t363 = s26 + s34
      t364 = s25 * t83
      t365 = s13 * (t61 + s25)
      t366 = t43 * t354
      t367 = t365 - t366
      t368 = t20 * t30
      t369 = t368 + s25
      t370 = t49 * s13
      t371 = t30 - t82
      t372 = s13 * (-t16 + t368) - t354 * t371
      t373 = s25 + s56
      t374 = t20 * t351
      t375 = t354 * (t16 + t374) - t342 * (t20 * t373 + s26)
      t376 = (t10 * t30 + t49) * s13 - t354 * (t368 - s25)
      t377 = t301 * s56
      t378 = t289 * s25
      t249 = t13 * (s46 * s56 + t304) - t20 * (s13 * t351 - s26 * s45 + 
     &(-s26 - s45 - s56) * s56) + s25 * (t249 + t377 + s26) + s26 * s46 
     &- t378
      t379 = intHL0s25s260s56x1021D2eps1()
      t380 = t9 * s26
      t381 = epinv * s56
      t382 = s12 + s13 - s16
      t383 = t259 - t247
      t384 = s16 - s34 - s13
      t385 = -t20 * t384 - t247 - t38
      t386 = t10 * t60 - t20 * t241 + s12 - t38
      t387 = s12 - s26
      t388 = t319 * s13
      t266 = (-s34 + t306) * s25 + (s25 - s12 - s13) * s26 + (-s25 - s12
     & + s26 - t342 + s45) * s45 + (-t387 - t266 + s46) * s46 + (-t98 - 
     &t387 + s56) * s56 + t388
      t306 = t382 * s25
      t389 = t343 * s34
      t390 = (-t20 * (s25 + s16 - s34 - s13 + s45) - s12 - t38 - s46) * 
     &s46
      t391 = t342 * s16
      t392 = intHs160000x0211D2eps1()
      t393 = t20 * (s13 - s45 - s46 - s56)
      t394 = -t393 + s25
      t395 = intHL0s25s260s56x1031D4eps1()
      t396 = intHs16s25s26s34s56x1112D2eps0()
      t397 = s26 - s56
      t398 = intHLs16s25s26s34s56x1311D4eps1()
      t399 = t43 + t100
      t400 = s26 * s34
      t263 = s13 * t263 - s16 * t399 + s25 * t43 + t10 * s34 * t83 + t20
     & * (s34 + s45 + s46) * s25 + t400 + t262 * s56
      t401 = intHL0s25s26s34s56x1220D4eps0()
      t402 = t260 + t43
      t403 = s25 + s26 + s34 + s56
      t404 = s13 * t402 - t399 * t403
      t405 = intHL0s25s26s34s56x1120D2eps1()
      t406 = t57 * s34
      t407 = epinv * s34
      t408 = t407 - s25 - s26 - s56
      t409 = -t20 * (s12 + s26 + s13 - s45 - s46 - s56) - t223
      t410 = t20 * t271
      t411 = s25 - s12 + s16 - s34
      t247 = t20 * t296 + s13 * t382 + (t247 - t342) * s26 + (-s12 + s16
     & - s34 + t332 + s45) * s45 + (t241 + t410 + s46) * s46 + (-t13 * t
     &60 + t411 + t69) * s56 - t276 + t334
      t276 = (t385 - s45) * s45
      t265 = -t20 * (s16 * t319 - (s12 + s13 + s16 + s26) * s34 - t296) 
     &- (t13 * t319 + s26 + t308) * s26 - (-t411 - t410 - s46) * s46 + (
     &-t10 * t265 - t13 * (s12 - s46) + t307) * s56 - t339 - t359 + t388
     & - (s16)**2 - t276
      t276 = -s26 * t383 - t343 * s34 - s56 * t386 - t276 + t297 - t306 
     &+ t388 - t390 - t391
      t319 = intHLs16s25s26s34s56x1231D6eps1()
      t332 = s25 + s45 + s46
      t334 = t20 * t387
      t339 = t10 * t332
      t288 = (-t288 - t49 - s26) * s26
      t257 = (-t257 + s25) * s25
      t382 = t20 * (t83 * (t325 + s34) + t305) + t297 + (t20 * t274 + s2
     &6) * s34 + (t13 * (s34 - s13) - s16 - t334 + t339) * s56 + t277 + 
     &t288 + t257
      t387 = t301 * s25
      t256 = t10 * (t363 + t82) * t83 + t20 * (s16 * t60 + t304 + t305) 
     &+ t297 - (t337 + t300) * s25 + (t10 * t274 + s26) * s34 + (-t13 * 
     &t256 + t215 + t295 - t334 + t387) * s56 - t324
      t411 = s16 + s26 + s56
      t412 = -t262 + t411
      t260 = -t260 + t411
      t244 = t244 - t100
      t413 = s13 * t412 + t244 * t260 + t49 * t260
      t414 = t9 * s34
      t415 = -t414 + s26 + s16 + s56
      t416 = t30 + t406
      t417 = t13 * epinv
      t418 = t417 - t10
      t419 = intHLs16s25s26s34s56x1321D6eps1()
      t420 = intHs160000x0121D2eps1()
      t421 = intHL0s25s26s34s56x1130D4eps1()
      t422 = intHs16s25s26s34s56x1113D4eps0()
      t423 = intHLs160000x0211D2eps1()
      t424 = intHL0s25s260s56x1012D2eps1()
      t425 = intHL0s25s260s56x1013D4eps1()
      t426 = intHL0s25s260s56x1010D0eps0()
      t427 = intHL0s25s260s56x1011D2eps0()
      t428 = intHL0s25s260s56x1020D2eps0()
      t429 = intHL0s25s26s34s56x1110D2eps0()
      t430 = intHs160s26s34s56x1013D4eps1()
      t431 = intHs16s25s26s34s56x1211D2eps0()
      t432 = intHL0s25s260s56x1022D4eps1()
      t433 = intHs16s25s26s34s56x1112D2eps1()
      t434 = intHL0s25s26s34s56x1220D4eps1()
      t435 = intHL0s25s260s56x1010D0eps1()
      t436 = intHL0s25s260s56x1020D2eps1()
      t437 = intHL0s25s26s34s56x1120D2eps0()
      t438 = intHL0s25s26s34s56x1210D2eps0()
      t439 = intHL0s25s26s34s56x1310D4eps0()
      t440 = intHs160s26s34s56x1013D4eps0()
      t441 = s26 - s13 + s45 + s56
      t320 = (t300 + t320 - s25) * s25
      t285 = t20 * ((-s25 + s12 + s26 + t224) * s56 + t285 - t305)
      t442 = t10 * (t277 + t362)
      t443 = t253 - t16
      t444 = s26 + s16 - s34 - s13 + s45 + s46
      t445 = -t259 + t82 + s26
      t217 = t83 * t217
      t446 = t445 * s34
      t384 = (t13 * t384 + t20 * (s12 + s45 + s46) + s25 + t33) * s56
      t274 = (t13 * t274 + s16 + s26 + t49) * s26
      t447 = t20 * (t305 + t217 - t309)
      t448 = t13 * (-t277 + t296)
      t449 = -t20 * t444 - s25
      t345 = -t284 + t345 + t296
      t450 = intHLs16s25s26s34s56x1221D4eps1()
      t451 = t20 * t343
      t452 = intHLs16s25s26s34s56x1141D6eps0()
      t453 = (s26)**2
      t454 = s12 * t43
      t455 = t332 * s16
      t456 = intHLs16s25s26s34s56x1131D4eps1()
      t457 = intHLs160000x0121D2eps1()
      t458 = t360 * s16
      t459 = t20 * (t458 - s25)
      t460 = t458 * t36
      t461 = t15 * (t26 * (t13 * t252 + s26 + t213 - t259 + t82) - t460)
     & - t17 * (t13 * (s13 - s46 - s56) - s26 - t213 + t459)
      t462 = intHLs160000x0121D2eps0()
      t459 = t10 * t207 - t212 + t459
      t212 = t26 * (-t253 + t10 * (s45 + s56) + t212 + t82) - t460
      t463 = intHLs160000x0122D4eps0()
      t464 = intHs160s26s34s56x1011D2eps1()
      t465 = intHs160s26s34s56x1012D2eps0()
      t466 = intHs160s26s34s56x1021D2eps0()
      t467 = intHLs16s25s26s34s56x1231D6eps0()
      t468 = t13 * t223
      t469 = s25 * s26
      t470 = intHs160000x0211D2eps0()
      t471 = intHL0s25s260s56x1031D4eps0()
      t472 = intHLs16s25s26s34s56x1211D4eps1()
      t473 = t10 * t218
      t474 = t301 * s46
      t475 = ddreal(8) * s45
      t476 = ddreal(8) * t83
      t477 = intHLs16s25s26s34s56x1121D4eps1()
      t252 = epinv * t226 - t252 * t301 - s26 - t275 + t289 - t475
      t478 = ddreal(1) / t43
      t479 = (t156)**2
      t480 = t122 * t418
      t481 = epinv * intHs160000x0111D0eps1()
      t482 = t24 * t52
      t483 = ddreal(1) / (ecossin)**2 * propW34 * im
      t484 = (intHL0s25s260s56x1011D2eps1() + intHL0s25s26s34s56x1110D2e
     &ps1()) * t21
      t485 = t52 * t478
      t486 = t25 * t52
      t487 = (t465 + t466) * s34
      t90 = t9 * t90
      t124 = t483 * (t28 * (t478 * (t233 * (-s34 * t437 * t52 - t232 - t
     &235) - t379 * (-t20 * (s13 * (-t95 - s26 - s56) + t356) + s25 * (-
     &t95 - t20 * (s45 * t9 + s46 * t9) - t99 + t381 - t380)) * t52) + t
     &484 * t233 * t52 * t478 + t485 * (t233 * (t435 - t436) + t404 * t4
     &34) * t9 + t486 * (s25 + s26 + s16 - s34 - s13 - s46 + s56 + t49))
     & + t482 * (t228 * (s13 + s26 - s56) + s34 * t394 * (t108 * t433 + 
     &t396)) + t52 * (t18 * (t9 * (-t105 * t266 - t124 * t413) + t94 * (
     &epinv * t266 + s26 * t383 + (t385 - s45) * s45 + s56 * t386 - t388
     & + t306 + t389 + t390 - t297 + t391) + t487 * t246) + t24 * (t246 
     &* (t21 * t464 + t349) + (t157 + t154) * (t20 * (s12 + s13 - s45) +
     & s26 + s16 - s34 - s46 - s56))) * t156 - t90 * t482 * t265 * t479)
      t40 = -t102 * (t20 * (-(s25 + s26 - s16) * s45 - t277 - t296 + t30
     &9 + t310) + (t308 - t307) * s46 + (-t20 * (s25 - s12 - s16 - s13 +
     & s45) - s46) * s56 - t311) - t113 * (epinv * t225 + ddreal(8) * t2
     &07 - t275 - t474) + t146 * t313 - t2 * t225 - t226 * t251 - t58 * 
     &(t20 * ((t271 + s46) * s46 - t269) - (-s46 * t20 + s25) * s56 + t2
     &67) - t59 * (t10 * (t312 * t83 + t296 + t304) - t378 - epinv * t31
     &3 + (t16 - t342) * s26 + (ddreal(9) * s25 + t10 * (s26 - s13 + s45
     & + s46)) * s56) + t96 * (s46 * t280 + t286 + t287) - t450 * (t10 *
     & t328 * t83 + t297 + t370 + epinv * (-s34 * t327 + (-t329 + s56) *
     & s56 - t323 - t324 + t330 + t331) + (-t253 + t16 + s26) * s34 + (t
     &295 + t339 - t259 - t334 - s16) * s56 + t257 + t277 + t288) - t472
     & * (epinv * t227 - s25 + t40 + t468 + t473 - t474 - t475 + t69) - 
     &t477 * t252
      t40 = t483 * (t24 * (t145 * t431 + t63 * (epinv * t145 - s13 - s16
     & - s26 + s34 + s56 - t230)) + t52 * (t110 * (t15 * t212 - t17 * t4
     &59) + t112 * (t15 * (epinv * t212 + t460) + t239 * (-t417 * t459 +
     & t460)) + t18 * (t101 * t222 + t118 * (t442 + epinv * (t20 * (t305
     & + t362 - t453) - (t273 + s13) * s25 - (t20 * (s12 - s13) + s16 + 
     &t16) * s26 + (-t468 - t73 + s25) * s56 + t446 + t259 * s16) + t317
     & - t320 - t285) + t121 * (epinv * t255 + s16 + s26 - s34 - t16 + t
     &289 - t377 - t476 + t49) + t123 * (epinv * t250 + s16 - s34 - t16 
     &+ t334 + t473 - t476) + t145 * t236 + t250 * t51 + t255 * t7) + t2
     &8 * t40 + t461 * t462))
      t207 = s34 * t233
      t80 = t107 * (-t20 * t345 + epinv * (-s34 * t333 - t335 + t336) - 
     &s26 * t343 + s56 * t449 - t318 + t389) + t115 * t341 + t116 * (t13
     & * (s56 * t223 + t469) + t20 * (s12 * t80 + s26 * t351 - s34 * t33
     &2 - t364 + t455) + s13 * (-t468 - t69 + s25) - s25 * s56 + s26 * t
     &223) + t127 * (-t445 * s34 + t274 + t384 - t447 + t448) + t211 * t
     &3 + t216 * t220 + t54 * (t10 * (-t284 + t340) + t297 + t298 + epin
     &v * t341 + (-t290 + s26) * s26 - t302 * s34 + (-t303 + t294 + t293
     & + t82) * s56)
      t212 = t133 * (s34 * t333 + t335 - t336) + t136 * (t10 * (t277 + t
     &217 - t296) + t370 + epinv * (-t448 - t274 + t446 - t384 + t447) -
     & (t300 + s13 + t49 - s25) * s25 + s26 * t443 - t443 * s34 + (-t10 
     &* t444 + s25) * s56) - t414 * t243 * t246 * t156
      t217 = epinv * t457
      t225 = t483 * t52
      t80 = t225 * (t18 * t212 + t24 * t80 + t28 * (-t1 * t208 - t117 * 
     &t256 - t120 * (t442 + epinv * (-s46 * t280 - t286 - t287) + t317 -
     & t320 - t285) + t209 * t214 - t221 * t227 + t326 * (s34 * t327 + (
     &t329 - s56) * s56 + t323 + t324 - t330 - t331) + t42 * (epinv * (-
     &t20 * (s13 * (s25 - s16 + s56) + (s25 + s16 - s34 - s45) * s26 + (
     &-s25 + s16 - s45) * s56 - t269 - t296 + t453 + t454 + t455) + s25 
     &* s34 + s46 * t43) - t365 + t366) + t44 * (-epinv * (t20 * (s25 * 
     &t321 + s46 * t441 + t268) + t229 * s25) + s25 * t145) - t207 * t43
     &8 * t478) + t217 * t461)
      t212 = epinv * t392
      t226 = t9 * t432
      t227 = t9 * t156
      t60 = t24 * (-t134 * (-t20 * t449 * s56 + epinv * t291 + s26 * t45
     &1 - t451 * s34 + t10 * t345 + t298) - t394 * t470 + t227 * (t153 *
     & t266 + t247 * t47 + t276 * t48)) - epinv2 * (t15 * (-t13 * (t361 
     &+ t296) - t20 * (t83 * (t38 + t223 + t82) - t277) + (propZ25 * (-s
     &56 * t281 + t277 + t283 - t309 - t310 + t362) + s16 + t337) * s25 
     &+ (s13 * t36 + s26 + t56 + t72) * s26 - (t307 + t342) * s34 + (t36
     & * t60 - t16 + t223 + t56 - t69) * s56) + t17 * (s13 * t223 + (s12
     & + t342) * s26 + (t98 - s25 + s12 - s26 - s56) * s56 - t83 * t358 
     &+ t359 - t361)) * t156 * t478
      t72 = t19 * t24
      t60 = t483 * (t18 * (-t128 * t394 + t52 * (-t130 * t372 - t143 * t
     &263 - t150 * t299 - t231 * t6 + (t75 * t9 - t140) * (t10 * t364 + 
     &t20 * (s34 * t353 + t304) + s13 * (-t16 + t224) - s16 * t354 + t36
     &3 * s25 + t16 * s56))) + t28 * (t52 * (t138 * t375 - t139 * t264 -
     & t144 * t270 + t148 * t249 + t149 * t278 + t155 * ((t213 - t259 + 
     &t82) * s25 - epinv * t270 + (t275 - s26) * s46 + (t16 + s46) * s56
     & + t469) - t394 * t471 + t456 * (-epinv * (t20 * (-s45 * t43 - t40
     &0 + t454) + (s16 - s34 - s46) * s25 + (s25 + s16 - s46) * s26 + (-
     &s25 + s16 - s46) * s56 - t296 + t355 + t453) - t20 * t367) - t226 
     &* t357 * t478) - t234 * t233 * t478) + t52 * t60 + t72 * (epinv * 
     &t394 - s16 + s25 - s26 + s34 + t215 - t219 - t49))
      t213 = t9 * t430
      t229 = t9 * t238
      t207 = t225 * (t18 * (-t131 * t376 - t151 * t292 - s34 * t246 * t1
     &56 * (t213 - t440)) + t28 * (-t145 * t62 - t248 * t314 - t315 * (-
     &t325 * s45 - s46 * t307 + (-s25 + s16 + s13 - s46 + t322 - s56) * 
     &s56 - t311 + t323 + t324) + t382 * t467) + t207 * t28 * (t229 - t4
     &39) * t478 - t126 * t24 * t291)
      t60 = t207 + t60 + t483 * (t52 * (t28 * ((t338 * t346 + t398 * (ep
     &inv * t346 - s16 + s25 - s26 + s34 + t215 - t219 - t49)) * s34 + t
     &350 * t357 * t478) - s34 * (t141 * (epinv * t409 + s26 + t215 - t2
     &19 + t82) + t409 * t422) * t18) + t9 * (t28 * (t52 * (t145 * t421 
     &+ t248 * t419 - t249 * t81 + t264 * t71 - t278 * t84 - t319 * t382
     & - t375 * t68 + t394 * t395) + t233 * t425 * t478) + t206 * (t132 
     &* t231 + t263 * t77 + t292 * t88 + t299 * t86 + t372 * t64 + t376 
     &* t66)) + t482 * (-t156 * (t247 * t92 + t266 * t29 + t276 * t93) -
     & t212 * t394))
      t207 = t9 * intHLs160000x0122D4eps1()
      t215 = intHLs160000x0111D2eps1() * t21
      t219 = (-t207 + t463) * s16 + t347 + t215
      t230 = t9 * intHLs16s25s26s34s56x1411D6eps1()
      t246 = t9 * intHLs16s25s26s34s56x1141D6eps1()
      t231 = t483 * (t52 * (t18 * ((s13 * t369 - t30 * t354 + t370) * (-
     &t180 + t129) + (-t182 + t142) * t145 * s34) + t28 * (t367 * (-t246
     & + t452) + (t178 - t147) * (-s46 * t397 + t20 * s25 * (t353 - s13)
     & + t304) + (-t230 + t348) * t231 * s34)) + t18 * (-t183 + t125) * 
     &t394 + t360 * t52 * t219 * t240)
      t247 = t9 * intHLs160000x0113D4eps1() - intHLs160000x0113D4eps0()
      t248 = ddreal(2) / ddreal(3)
      t249 = ddreal(1) / ddreal(3)
      t250 = epinv * t30
      t252 = s25 * t108
      t255 = -s16 + s25 - s26 + s34 - s56
      t257 = epinv * s25
      t259 = t9 * t48
      t263 = t9 * t425
      t264 = epinv * t456
      t267 = t117 * t403
      t268 = t9 * t477
      t269 = t424 * t9 * t478
      t270 = t420 * epinv
      t271 = t349 * t156
      t274 = t113 * t418
      t276 = s34 * t396
      t277 = -t435 + t436
      t278 = t379 * (t95 + s26 + s56)
      t280 = t94 * (t257 + s26 + s16 - s34 + s56)
      t281 = s34 * t433
      t283 = t118 * t24
      t284 = t18 * t107
      t285 = t108 * (t281 * t18 + t283 * t210) + t156 * (t18 * (-t255 * 
     &t92 + t280 + t487 - t90) + t24 * ((-t173 + t5) * s34 + t21 * (t91 
     &+ t464))) + t28 * ((-t108 * t35 + t478 * (-t401 + t437 + t438)) * 
     &s34 + t478 * (t277 * t9 - t350 * t352 - t278) - t484 * t478) + (-t
     &28 * t42 + t284) * t258 * epinv
      t258 = t52 * t285 + t18 * (t156 * (t52 * (t9 * (-s34 * t243 + t255
     & * t47) + t407 * t4) + t103 + t104) + t276 * t52) + t28 * (t478 * 
     &(t232 + t235) + t52 * (-t102 * t258 + t120 * t210 - t25 + t274 - t
     &462) + t146 + t485 * ((t434 - t237) * s34 + t352 * t432) * t9 + t5
     &2 * (s34 * t405 * t478 - t457) * epinv) + t53 * (t93 + t470) + t28
     & * (t52 * (-t30 * t326 - t450 * (t250 - s34) + t478 * (t426 - t427
     & - t428 - t429) + t59 * (t252 + s26 + s56) + t267 - t268) - t269) 
     &+ t52 * (t15 * (t12 * t216 + t22 * (t228 + t7) + t251 * t26) + t18
     & * (t101 * t260 + t133 * t258 + t270)) + t482 * (-t121 * t9 - t127
     & * t371 - t136 * (epinv * t371 + s16 + s26 - s34 + s56) + t54 * (-
     &epinv * t261 + s34) + t480 + t271) + t482 * (-t115 * t261 + t116 *
     & t210)
      t261 = t18 * t152
      t154 = t24 * t154
      t285 = t105 * t9
      t286 = t285 * t18
      t177 = t177 * t15
      t287 = t162 * t403
      t288 = epinv * t111
      t224 = t156 * (t52 * (t261 + t154 + epinv2 * (t15 * ((-propZ25 * t
     &210 + t10) * s25 + t351 * t36 + t224) + t17 * (t241 + t374)) * t47
     &8) - t286) + t52 * (-t159 * t18 * t260 + t177) + t258 + (-t172 + t
     &288 - t287) * t52 * t28
      t241 = t52 * t15
      t258 = t85 * (t251 + t216 - t228 + t1 - t7) + t241 * (t22 * t8 + t
     &196)
      t260 = t183 - t125
      t289 = t260 * t24
      t72 = t13 * t258 + t20 * (t24 * (t156 * t29 + t128 + (t108 * t141 
     &+ t156 * t440) * t52 * s34) + t28 * (-t52 * (s34 * t439 * t478 + t
     &110) + t315) + t9 * (t24 * (-t153 * t156 + t52 * (-(t156 * t430 + 
     &t77) * s34 - t210 * t75 - t255 * t64)) - t206 * t66 * t30 + t52 * 
     &((t238 * t478 - t71) * s34 - t210 * t84 - t352 * t81 - t43 * t68) 
     &* t28) + t28 * (t155 + t264) + t52 * (t18 * (t131 * t30 + t212 - t
     &259) + t24 * ((t143 + t422) * s34 + t140 * t210 + t255 * (t126 + t
     &130)) + t28 * ((t398 + t139) * s34 + t138 * t43 + t148 * t352 + t1
     &49 * t210 + t395 * t9 - t471)) + t28 * (t263 - t234) * t478 + t72 
     &* t108 + t52 * (t134 * t24 * t255 - t112 * t28) * epinv) - t36 * (
     &-t186 + t289) + t224 + ddreal(9) * t85 * (t2 - t8)
      t212 = s25 + s26 - s34 + s56
      t224 = t36 * propW16
      t238 = t308 * propW16
      t258 = t224 * s16
      t290 = t15 * (t11 + t258 + t20) + t17 * (ddreal(1) + t238)
      t291 = (t392 + t420) * epinv
      t292 = epinv2 * (t15 * (t224 * t30 - t10 + t11) + t17 * (t368 * pr
     &opW16 - ddreal(1))) * t156
      t293 = t110 * t290
      t294 = -t5 + t173
      t295 = -t47 * t9 + t92
      t243 = t9 * t243
      t255 = t295 * t255
      t91 = (-t91 - t464) * t21
      t290 = t112 * (epinv * t290 - t258 * t240)
      t296 = t42 * t43
      t212 = t156 * (t52 * (-t261 - t154) + t286) + t156 * (t18 * (t52 *
     & ((-epinv * t4 + t243 - t465 - t466) * s34 - t280 + t90 + t255) - 
     &t103 - t104) + t482 * (s34 * t294 + t91)) + t52 * (-t284 * epinv *
     & t212 - t283 * t43 * t108 + t290) - t53 * t93 + t52 * (t15 * (t12 
     &* (-t216 - t470) + t22 * (-t7 - t8)) + t18 * (-t133 * t212 - t281 
     &- t291) + t24 * (t169 * t9 + t43 * (-t115 - t116) - t54 * (epinv *
     & t43 + s34) - t271) + t28 * (s34 * t58 - t120 * (epinv * t411 + s2
     &5 + s26 + s56) + t146 * t344 - t411 * t96 - t59 * (t407 + t95 + s2
     &6 + s56) - t267) + t292 + t293) + t52 * (epinv * ((t35 + t44) * s3
     &4 + t296) + t25) * t28 - t177 * t52 + (t102 * t43 - t171 + t172 + 
     &t287 - t288) * t52 * t28
      t213 = t213 - t440
      t267 = t71 * t9 - t139
      t283 = t9 * t81
      t284 = t283 - t148
      t286 = t9 * t84
      t287 = t286 - t149
      t297 = t9 * t68
      t298 = (-t29 + t200) * t156
      t299 = propW16 * s16
      t302 = t299 * t240
      t304 = t302 * t52 * t247
      t305 = t20 + epinv
      t306 = t20 * t441 + t16
      t309 = t10 * s56
      t310 = t309 + t16 + t410
      t311 = -t20 * (s16 - s34 + s13 - s45 - s56) + s25 - t56
      t301 = t272 * t36 + t279 * t301 + s25 + t99
      t262 = t20 * (s12 + s26 + s16) - t262
      t282 = -t282 + s25
      t312 = s12 + s26 + s16 - s34 + s56
      t313 = t10 * t312
      t312 = t20 * t312
      t317 = -t312 + s25
      t273 = -t316 + t38 + t373 + t273
      t318 = -s13 + s25 + s26 + s45 + s56
      t320 = t10 * t272
      t321 = t37 * s26
      t322 = t299 * t10
      t299 = ddreal(12) * t299
      t323 = t15 * (t11 + t299 + t20) + t17 * (ddreal(1) + t322)
      t66 = t66 * t9
      t324 = t66 - t131
      t86 = t9 * t86
      t325 = t9 * t132
      t327 = s34 * t338
      t193 = (t193 + t194) * t9
      t194 = t19 * t305
      t88 = t9 * t88
      t328 = t52 * (t43 * (t102 - t264) + t172 - t201 - t288 - t74) * t2
     &8
      t91 = t18 * (t52 * (t88 * t273 - t175) + t227 * t105) - t52 * (t15
     &4 * t156 + t177) + t156 * (t18 * (t52 * ((t243 - t465 - t466) * s3
     &4 + t255 - t280 + t90) - t103 - t104) + t91 * t482) - t24 * (t194 
     &+ t126 + t127 + t128) + t52 * (t18 * ((-t392 - t420) * epinv - t14
     &3 * t262 - t151 * t273 + (t86 - t150) * (t279 * t37 + t320 + t82 +
     & t99) + t324 * (t313 + s25) - t281 + t325) + t24 * (-t130 * t317 -
     & t133 * t272 - t134 * ((-t20 + epinv) * s25 - t313) - t136 * (t257
     & - t312) - t140 * t282 - t271) - t28 * (t155 * (t57 * t318 + s25) 
     &+ t315 * t43 + t327) + t290 + t292 + t293 - t23 * t51) + (t24 * (t
     &282 * t75 + t317 * t64) + t18 * t77 * t262) * t52 * t9 + t52 * (ep
     &inv * ((-t398 + t35) * s34 + t296) + t193) * t28 + t328
      t154 = epinv * t423
      t175 = t403 * t450
      t34 = (-t162 + t117) * t34
      t177 = t70 * (t115 + t116)
      t169 = t18 * (t52 * (t156 * (-t294 * t412 + t4 * t416) - t159 * t3
     &01) + t63) + t53 * (t6 + t216 + t470) + t52 * ((t106 * t18 - t323 
     &* t457) * epinv + t15 * (t12 * t3 + t26 * (t221 + t251)) + t24 * (
     &t107 * (t20 * (t108 * t272 + s26 - s34 + s56) + s25) + t118 * (-t1
     &0 * (s12 + s16 - s34) - t20 * t380 + t252) + t169 * t418 + t54 * (
     &-t20 * (-epinv * s26 + s25) + t257 + t303 - t320 - t99 - t321) + t
     &176) + t28 * (t120 * (epinv * t311 + s25 + t69) + t146 * t310 - t3
     &06 * t58 + t311 * t96 + t326 * t403 - t472 * t9 + t59 * (epinv * t
     &310 + s25 + t69) - t154) + t101 * t18 * t301) + t177 * t482 + t52 
     &* ((-t306 * t44 + t175) * epinv - t209 - t25 - t268 + t274 + t34) 
     &* t28 - t462 * t323 * t52
      t243 = -t259 + t93
      t90 = t90 - t152
      t252 = (t120 - t44 + t59) * epinv + t146 - t58 + t96
      t255 = epinv * t155 + t144
      t262 = s25 - s13 + s45 + s46 + t61
      t268 = t108 * s12
      t271 = t20 * (t268 + s25 + s26)
      t272 = s25 - s26 + s56
      t273 = s16 + s56
      t274 = -t98 + t99 + t307
      t140 = t140 + t130 + t19
      t280 = t9 * t319
      t281 = t398 * s34
      t282 = t156 * t52
      t27 = t18 * (t9 * (t75 + t64) - t203 * t52 + t52 * (t86 - t150) * 
     &t87 + t52 * (t88 - t151) * t89) + t24 * (t52 * (t325 - t199 - t198
     & - t163) - t126 + t298 - t282 * t243 * t50 - t282 * t295 * t31) + 
     &t28 * (t486 - t315) + t52 * (t28 * (-t144 * t274 - t146 * t262 - t
     &155 * (epinv * t274 - s25 + s26 - s56) - t210 * t267 + t272 * (-t9
     & * (t81 + t84) + t148 + t149) + t421 * t9 - t456 * (t55 - t33) + (
     &-t280 + t467) * (t321 + t16 + t273 + t49) + (-t419 * t9 + t314) * 
     &(t38 + t32 + s16 + s34 + s56) + t221 + t281) - t27 * t62 + t210 * 
     &(t77 * t9 - t143) * t18) + t23 * t140 - t241 * (t3 + t6) * t22
      t33 = t52 * ((t161 - t423) * epinv + t114 * t472 + t21 * t477 - t3
     &26 * t50 + t46 * t58 + t9 * (t113 + t434) + t110 + t251) * t28
      t55 = t15 * (t22 * (-t4 - t5) + t26 * (-t1 - t2)) + t24 * (-t115 *
     & t70 + t173 - t174) + t28 * ((-t111 + t112) * epinv + t172) - t261
     & * t156
      t12 = t52 * t55 + t24 * (-t108 * t63 + t52 * (-t156 * t158 - t164 
     &- t165 - t270) - t431) - t53 * t7 + t52 * (t15 * (-t12 * t236 - t2
     &16 * t22 - t26 * (t209 + t401 + t405)) + t18 * (-t114 * t121 - t11
     &6 * t242 - t118 * (epinv * t358 + t271) - t123 * t418 - t127 * (t3
     &0 + t49) - t136 * (t250 + t16 + t109) + (t159 - t101) * (t16 + s26
     & - s16 + s34 - s56 + t49) - t481) + t28 * (t102 * t272 + t120 * (e
     &pinv * (t98 - t99 + t279) + t271) - t42 * (s26 * t305 - t381 + t95
     &) + t44 * (epinv * t46 - s25) - t450 * (-t49 * t9 - s16 - s56 - t1
     &6 + t257 - t321) - t59 * (t20 * (epinv * t262 - s26) - t16) + t96 
     &* (-t99 + t279 + t245) + (-t162 + t117) * (t275 + t38 + s56 + t49)
     &)) - t241 * t228 * t22 - t166 * t156 + t33
      t12 = -t10 * t52 * (t31 * t18 * t324 + (-t297 + t138) * t307 * t28
     &) - t13 * (t241 * t51 * t22 + t85 * (-t401 - t405 - t236 - t209 + 
     &t216 + t228 - t1 - t2 + t4 + t5)) - t20 * t27 + t36 * (t52 * (t18 
     &* (t184 + t185) + t28 * ((t230 - t348) * s34 + t179 * t397 + t70 *
     & (t246 - t452))) + t85 * (s25 * t140 - t3 - t6 + t62)) - t12 + ddr
     &eal(9) * t85 * t51
      t26 = ddreal(1) + t57
      t27 = t15 * (t11 - t299 + t20) + t17 * (ddreal(1) - t322)
      t31 = t15 * (t11 - t258 + t20) + t17 * (ddreal(1) - t238)
      t33 = t308 + t43
      t38 = t20 * (s25 - s12 - s16 + s34) + s56 - t38
      t46 = t253 - t99 - t100 - s26
      t53 = t10 * t279 + t20 * t273 + s25
      t55 = t20 * t43
      t75 = t275 + t99 - s26
      t39 = t39 + t309 + t69
      t77 = s25 * t36 + t10 * (s34 + s56) - t300 - t49
      t86 = t20 * (s26 + s45 + s46 + s56) + s25 - t253
      t87 = -t20 * t254 + s25
      t73 = t316 - t300 - t73 + s25
      t88 = t337 - t99 + t97 - t100
      t89 = -s13 + s25 + s45 + s46 + s56
      t97 = -s12 - s13 + s25 + s45 + s46 + s56
      t98 = s12 - s13 + s25 + s45 + s46 + s56
      t30 = t483 * (t28 * (t52 * (t114 * (t472 + t477) + t138 * t39 + t1
     &39 * t73 + t144 * t46 + t146 * (t20 * t218 - t83) + t148 * t75 + t
     &149 * t38 + t155 * (epinv * t46 - s26 + s56 + t16) - t154 + t2 - t
     &209 + t25) - t269) + t482 * (-t126 * t369 - t127 * t30 - t133 * (s
     &12 - s16 - s26 + s34 - s56) - t291 + t325) - t206 * (t205 + t204))
      t11 = t225 * (t18 * (-t136 * (-t250 + t32) - t76 - t79) + t24 * (-
     &t134 * (epinv * t369 - t56 - t82) - t276) + t478 * (epinv2 * (t15 
     &* (t224 * t43 + t11 + t20) + t17 * (t55 * propW16 + ddreal(1))) + 
     &t28 * (t426 - t427 - t428 - t429)) - t217 * t31)
      t11 = t11 + t30 + t483 * (t28 * (t478 * (t232 + t235) + t52 * (t9 
     &* (-s34 * t237 * t478 - t38 * t84 - t39 * t68 - t75 * t81) + t405 
     &* t408 * t478) + t485 * (t434 * t9 - t401) * t403 + t485 * (t226 -
     & t350) * t352) + t52 * (-t18 * (t197 + t236 + t202 + t481) - t31 *
     & t462) + t482 * (t227 * (t153 * t97 + t48 * t98) - t6 - t216 - t47
     &0) - t194 * t18) + t483 * (-t128 * t18 + t24 * (t282 * (-t29 * t97
     & - t89 * t92 - t93 * t98) - t431) + t485 * (-t484 - t278) * t28 + 
     &t52 * (-t433 * epinv * t24 + t28 * t478 * (t437 + t438) - t141 * t
     &26 * t18) * s34 + t52 * (t18 * (t192 + t187 + t188 + t189 + t190 +
     & t191) + t28 * (t277 * t478 - t71 * t73) + t47 * t24 * t89 * t156)
     & * t9)
      t17 = t225 * (t18 * (-t101 * t41 - t116 * t87 - t118 * (epinv * t8
     &7 + t16 + t56 + t69)) + t24 * (-t107 * (-t20 * (-t268 + t250) + s2
     &5) - t115 * t53 - t54 * (epinv * t53 - s26 - s34 + s56 - t56 - t82
     &) + t173 - t174) + t28 * (t102 * t33 + t117 * t77 + t120 * (epinv 
     &* t88 - t10 * t223 + t16 - t69) + t42 * (epinv * t33 - t55) + t44 
     &* (epinv * t86 - t82) + t58 * t86 + t59 * (t10 * (epinv * s13 + s5
     &6) - t20 * (epinv * (t83 + t61) - s26) + t387) + t96 * t88 + t171 
     &+ t172 - t288))
      t30 = t483 * (t28 * (t478 * ((-t229 + t439) * t52 * s34 + t234 - t
     &263) + t52 * (-t314 * t402 + t9 * (t402 * t419 - t395 - t421) - t1
     &75 - t221 - t251 + t471 + t62) + t52 * (t280 - t467) * (t55 + s34)
     &) + t206 * (s34 * t422 + t51))
      t11 = -t10 * t483 * (t18 * (t52 * (t182 - t142) * s34 - t125 + t18
     &3 + t52 * (t180 - t129) * t50) + t52 * (-t178 + t147) * t272 * t28
     &) - ddreal(4) / ddreal(3) * t11 + t20 * t225 * (t170 * t18 + t24 *
     & (t163 + t3)) - t248 * (t17 + t483 * (t24 * (-t26 * t63 - t52 * (t
     &157 * t156 + t4 + t5)) + t52 * (-t110 * t27 - t112 * (epinv * t27 
     &+ ddreal(12) * t302) + t28 * (t161 * t21 - t162 * t77 - t1)) + t20
     &6 * (-t123 * (epinv * t10 - t37) + t418 * (-t121 - t122) + t479 * 
     &t90 * (s12 + s13 + s16 + s26 - s34 - s45 - s46) + (t285 - t103 - t
     &104) * t156 * t97 + t94 * (t89 * t9 + t268) * t156 + t159 * t41)))
     & + ddreal(8) / ddreal(3) * t30 - ddreal(16) * t302 * t225 * t247 -
     & ddreal(16) / ddreal(3) * t225 * t28 * (t43 * t456 + t281) + ddrea
     &l(8) * t225 * (t28 * ((t230 - t348) * s34 + t43 * (t246 - t452)) +
     & t219 * t240 * propW16)
      t17 = epinv * t70
      t26 = t56 + t393 + s25
      t27 = t393 + s25
      t6 = -t51 + t6
      t16 = t18 * (-t121 * t21 - t173 - t7 - t8) + t28 * ((-t423 - t457)
     & * epinv - t113 * t21 + t9 * (-t472 - t477) - t2 - t209 + t221 + t
     &25 + t251 - t462) + t18 * (-t118 * (t56 + t16 + t69 - t17) - t122 
     &* t21 - t54 * (t56 - t17 + s26 + s34 - s56 + t82) + t164 - t167 - 
     &t176 + t177) + t28 * (-t120 * (-epinv * t26 + s25 + t69) + t146 * 
     &t27 + t26 * t96 + t403 * (epinv * t450 + t326) - t59 * (-epinv * t
     &27 + s25 + t69) - t160 - t34) - t163 * t24 + t23 * (-t3 - t4 - t5)
     & - t63 * s25 * t24
      t6 = struc60PM * (-t10 * t28 * t255 - t20 * (t252 * t28 + (t9 * (t
     &153 + t47 + t48) - t29 - t92 - t93) * t156 * t18) - t156 * (t24 * 
     &(t9 * (t94 + t105) - t103 - t104) + t90 * t156 * t18)) + struc9PM 
     &* (-t10 * t28 * (-t144 * t45 + t307 * (t283 - t148) + t78 * (t286 
     &- t149)) - t13 * t239 * t168 - t20 * (t18 * (t134 * t135 - t203 - 
     &t204 - t205 - t76 - t79) + t28 * (t327 - t201 - t74) + t18 * ((t19
     &5 + t126 + t127 + t128) * s25 + t136 * t137 + t65 * (t64 * t9 - t1
     &30) + t67 * (t66 - t131) + t9 * (-t132 + t189 + t190 + t191 + t192
     &)) + t28 * (-t155 * (t57 * t45 + s25) + t43 * (t315 + t264) + t193
     & + t407 * t398) - t23 * t6 - t198 * t24) - t36 * (t18 * (s25 * t26
     &0 - t184 - t185) + t239 * t6 - t186 * s25) - t16)
      t15 = t483 * (t72 * struc2PM + (struc3PM + struc8PM) * (t13 * t85 
     &* (t1 - t7 - t8 + t216 + t470) - t20 * (t24 * (t298 - t126 - t127 
     &- t128) + t52 * (t28 * (t267 * s34 + t144 * t344 - t155 * (t407 + 
     &t95) + t210 * t287 + t352 * t284 + t43 * (t297 - t138)) - t196 * t
     &15 + t259 * t18 + s34 * (t156 * t213 - t141) * t24) - (t134 + t136
     & + t19) * t24 * epinv) + t36 * (t85 * t2 + t186 + ((-t207 - t462 +
     & t463 - t217) * s16 + t215 + t347) * t52 * t240 * propW16) - ddrea
     &l(12) * t304 - t212) + (-t20 * t255 - t252) * t28 * struc61PM * t5
     &2)
      t4 = (-t10 * t225 * t458 * t240 * t247 + t20 * t231 + t249 * (t80 
     &+ t40 + t124 + t483 * (t28 * (t478 * (t52 * (t233 * (t414 * t237 -
     & t426 + t427 + t428 + t429) - t401 * t404 - t405 * (s13 * (t43 - t
     &406) + t399 * t408)) + t9 * t233 * t424) + t52 * ((t119 * t208 + t
     &35 * (epinv * t208 + s25 - t393)) * s34 + t162 * t256) + t52 * (-t
     &111 * t208 + t214 * t423) * epinv) + t206 * (t156 * (t266 * (t103 
     &+ t104) + t4 * (s13 * t416 + t244 * t415 + t49 * t415) + t413 * t5
     &) + t480 * t394 + t481 * t145 - t159 * t222) + t482 * ((t106 * t21
     &1 + t220 * t420) * epinv + t152 * t265 * t479))) - t248 * t60 + t4
     &83 * t206 * t8 * t394) * struc1PM + t11 * struc7PM
      result = t248 * t483 * (t14 * struc10PM + t52 * t6 + (t10 * (t52 *
     & (t18 * t243 + t28 * (t144 * t318 - t284 * t307 - t287 * t78)) + (
     &-t213 * t52 * s34 - t200 + t29) * t156 * t24) + t13 * (t85 * (t251
     & + t221 + t3) + t241 * (t22 * (t7 + t8) + t196)) - t20 * t91 - ddr
     &eal(24) * t304 - t36 * (-t186 + t289 + t482 * (-t141 - t142 + t182
     &) * s34 + t85 * (t51 - t1) - t206 * t181 * (t312 + s25)) + t169 + 
     &ddreal(12) * t52 * t240 * propW16 * t219 + ddreal(9) * t85 * (t2 -
     & t7 - t8)) * struc5PM + t12 * struc6PM) + ddreal(4) / ddreal(3) * 
     &t15 + t4

           ampNonresonantLightFullImC4PM = result
       end function ampNonresonantLightFullImC4PM

       function ampNonresonantLightFullImC4PP()
           implicit none
           type(dd_complex) :: ampNonresonantLightFullImC4PP

           type(dd_complex) :: result

      result = ddreal(0)

           ampNonresonantLightFullImC4PP = result
       end function ampNonresonantLightFullImC4PP

       function ampNonresonantLightFullImC7MM()
           implicit none
           type(dd_complex) :: ampNonresonantLightFullImC7MM

           type(dd_complex) :: result

      result = ddreal(0)

           ampNonresonantLightFullImC7MM = result
       end function ampNonresonantLightFullImC7MM

       function ampNonresonantLightFullImC7MP()
           implicit none
           type(dd_complex) :: ampNonresonantLightFullImC7MP

           type(dd_complex) :: result

      result = ddreal(0)

           ampNonresonantLightFullImC7MP = result
       end function ampNonresonantLightFullImC7MP

       function ampNonresonantLightFullImC7PM()
           implicit none
           type(dd_complex) :: ampNonresonantLightFullImC7PM

           type(dd_complex) :: result

      result = ddreal(0)

           ampNonresonantLightFullImC7PM = result
       end function ampNonresonantLightFullImC7PM

       function ampNonresonantLightFullImC7PP()
           implicit none
           type(dd_complex) :: ampNonresonantLightFullImC7PP

           type(dd_complex) :: result

      result = ddreal(0)

           ampNonresonantLightFullImC7PP = result
       end function ampNonresonantLightFullImC7PP

       function ampNonresonantLightFullMM()
           implicit none
           type(dd_complex) :: ampNonresonantLightFullMM
           type(dd_complex) ::  t1,t10,t100,t1000,t1001,t1002,t1003,t1004,t1005,t1006,t1007,t1008
           type(dd_complex) ::  t1009,t101,t1010,t1011,t1012,t1013,t1014,t1015,t1016,t1017,t1018,t1019
           type(dd_complex) ::  t102,t1020,t1021,t1022,t1023,t1024,t1025,t1026,t1027,t1028,t1029,t103
           type(dd_complex) ::  t1030,t1031,t1032,t1033,t1034,t1035,t1036,t1037,t1038,t1039,t104,t1040
           type(dd_complex) ::  t1041,t1042,t1043,t1044,t1045,t1046,t1047,t1048,t1049,t105,t1050,t1051
           type(dd_complex) ::  t1052,t1053,t1054,t1055,t1056,t1057,t1058,t1059,t106,t1060,t1061,t1062
           type(dd_complex) ::  t1063,t1064,t1065,t1066,t1067,t1068,t1069,t107,t1070,t1071,t1072,t1073
           type(dd_complex) ::  t1074,t1075,t1076,t1077,t1078,t1079,t108,t1080,t1081,t1082,t1083,t1084
           type(dd_complex) ::  t1085,t1086,t1087,t1088,t1089,t109,t1090,t1091,t1092,t1093,t1094,t1095
           type(dd_complex) ::  t1096,t1097,t1098,t1099,t11,t110,t1100,t1101,t1102,t1103,t1104,t1105
           type(dd_complex) ::  t1106,t1107,t1108,t1109,t111,t1110,t1111,t1112,t1113,t1114,t1115,t1116
           type(dd_complex) ::  t1117,t1118,t1119,t112,t1120,t1121,t1122,t1123,t1124,t1125,t1126,t1127
           type(dd_complex) ::  t1128,t1129,t113,t1130,t1131,t1132,t1133,t1134,t1135,t1136,t1137,t1138
           type(dd_complex) ::  t1139,t114,t1140,t1141,t1142,t1143,t1144,t1145,t1146,t1147,t1148,t1149
           type(dd_complex) ::  t115,t1150,t1151,t1152,t1153,t1154,t1155,t1156,t1157,t1158,t1159,t116
           type(dd_complex) ::  t1160,t1161,t1162,t1163,t1164,t1165,t1166,t1167,t1168,t1169,t117,t1170
           type(dd_complex) ::  t1171,t118,t119,t12,t120,t121,t122,t123,t124,t125,t126,t127
           type(dd_complex) ::  t128,t129,t13,t130,t131,t132,t133,t134,t135,t136,t137,t138
           type(dd_complex) ::  t139,t14,t140,t141,t142,t143,t144,t145,t146,t147,t148,t149
           type(dd_complex) ::  t15,t150,t151,t152,t153,t154,t155,t156,t157,t158,t159,t16
           type(dd_complex) ::  t160,t161,t162,t163,t164,t165,t166,t167,t168,t169,t17,t170
           type(dd_complex) ::  t171,t172,t173,t174,t175,t176,t177,t178,t179,t18,t180,t181
           type(dd_complex) ::  t182,t183,t184,t185,t186,t187,t188,t189,t19,t190,t191,t192
           type(dd_complex) ::  t193,t194,t195,t196,t197,t198,t199,t2,t20,t200,t201,t202
           type(dd_complex) ::  t203,t204,t205,t206,t207,t208,t209,t21,t210,t211,t212,t213
           type(dd_complex) ::  t214,t215,t216,t217,t218,t219,t22,t220,t221,t222,t223,t224
           type(dd_complex) ::  t225,t226,t227,t228,t229,t23,t230,t231,t232,t233,t234,t235
           type(dd_complex) ::  t236,t237,t238,t239,t24,t240,t241,t242,t243,t244,t245,t246
           type(dd_complex) ::  t247,t248,t249,t25,t250,t251,t252,t253,t254,t255,t256,t257
           type(dd_complex) ::  t258,t259,t26,t260,t261,t262,t263,t264,t265,t266,t267,t268
           type(dd_complex) ::  t269,t27,t270,t271,t272,t273,t274,t275,t276,t277,t278,t279
           type(dd_complex) ::  t28,t280,t281,t282,t283,t284,t285,t286,t287,t288,t289,t29
           type(dd_complex) ::  t290,t291,t292,t293,t294,t295,t296,t297,t298,t299,t3,t30
           type(dd_complex) ::  t300,t301,t302,t303,t304,t305,t306,t307,t308,t309,t31,t310
           type(dd_complex) ::  t311,t312,t313,t314,t315,t316,t317,t318,t319,t32,t320,t321
           type(dd_complex) ::  t322,t323,t324,t325,t326,t327,t328,t329,t33,t330,t331,t332
           type(dd_complex) ::  t333,t334,t335,t336,t337,t338,t339,t34,t340,t341,t342,t343
           type(dd_complex) ::  t344,t345,t346,t347,t348,t349,t35,t350,t351,t352,t353,t354
           type(dd_complex) ::  t355,t356,t357,t358,t359,t36,t360,t361,t362,t363,t364,t365
           type(dd_complex) ::  t366,t367,t368,t369,t37,t370,t371,t372,t373,t374,t375,t376
           type(dd_complex) ::  t377,t378,t379,t38,t380,t381,t382,t383,t384,t385,t386,t387
           type(dd_complex) ::  t388,t389,t39,t390,t391,t392,t393,t394,t395,t396,t397,t398
           type(dd_complex) ::  t399,t4,t40,t400,t401,t402,t403,t404,t405,t406,t407,t408
           type(dd_complex) ::  t409,t41,t410,t411,t412,t413,t414,t415,t416,t417,t418,t419
           type(dd_complex) ::  t42,t420,t421,t422,t423,t424,t425,t426,t427,t428,t429,t43
           type(dd_complex) ::  t430,t431,t432,t433,t434,t435,t436,t437,t438,t439,t44,t440
           type(dd_complex) ::  t441,t442,t443,t444,t445,t446,t447,t448,t449,t45,t450,t451
           type(dd_complex) ::  t452,t453,t454,t455,t456,t457,t458,t459,t46,t460,t461,t462
           type(dd_complex) ::  t463,t464,t465,t466,t467,t468,t469,t47,t470,t471,t472,t473
           type(dd_complex) ::  t474,t475,t476,t477,t478,t479,t48,t480,t481,t482,t483,t484
           type(dd_complex) ::  t485,t486,t487,t488,t489,t49,t490,t491,t492,t493,t494,t495
           type(dd_complex) ::  t496,t497,t498,t499,t5,t50,t500,t501,t502,t503,t504,t505
           type(dd_complex) ::  t506,t507,t508,t509,t51,t510,t511,t512,t513,t514,t515,t516
           type(dd_complex) ::  t517,t518,t519,t52,t520,t521,t522,t523,t524,t525,t526,t527
           type(dd_complex) ::  t528,t529,t53,t530,t531,t532,t533,t534,t535,t536,t537,t538
           type(dd_complex) ::  t539,t54,t540,t541,t542,t543,t544,t545,t546,t547,t548,t549
           type(dd_complex) ::  t55,t550,t551,t552,t553,t554,t555,t556,t557,t558,t559,t56
           type(dd_complex) ::  t560,t561,t562,t563,t564,t565,t566,t567,t568,t569,t57,t570
           type(dd_complex) ::  t571,t572,t573,t574,t575,t576,t577,t578,t579,t58,t580,t581
           type(dd_complex) ::  t582,t583,t584,t585,t586,t587,t588,t589,t59,t590,t591,t592
           type(dd_complex) ::  t593,t594,t595,t596,t597,t598,t599,t6,t60,t600,t601,t602
           type(dd_complex) ::  t603,t604,t605,t606,t607,t608,t609,t61,t610,t611,t612,t613
           type(dd_complex) ::  t614,t615,t616,t617,t618,t619,t62,t620,t621,t622,t623,t624
           type(dd_complex) ::  t625,t626,t627,t628,t629,t63,t630,t631,t632,t633,t634,t635
           type(dd_complex) ::  t636,t637,t638,t639,t64,t640,t641,t642,t643,t644,t645,t646
           type(dd_complex) ::  t647,t648,t649,t65,t650,t651,t652,t653,t654,t655,t656,t657
           type(dd_complex) ::  t658,t659,t66,t660,t661,t662,t663,t664,t665,t666,t667,t668
           type(dd_complex) ::  t669,t67,t670,t671,t672,t673,t674,t675,t676,t677,t678,t679
           type(dd_complex) ::  t68,t680,t681,t682,t683,t684,t685,t686,t687,t688,t689,t69
           type(dd_complex) ::  t690,t691,t692,t693,t694,t695,t696,t697,t698,t699,t7,t70
           type(dd_complex) ::  t700,t701,t702,t703,t704,t705,t706,t707,t708,t709,t71,t710
           type(dd_complex) ::  t711,t712,t713,t714,t715,t716,t717,t718,t719,t72,t720,t721
           type(dd_complex) ::  t722,t723,t724,t725,t726,t727,t728,t729,t73,t730,t731,t732
           type(dd_complex) ::  t733,t734,t735,t736,t737,t738,t739,t74,t740,t741,t742,t743
           type(dd_complex) ::  t744,t745,t746,t747,t748,t749,t75,t750,t751,t752,t753,t754
           type(dd_complex) ::  t755,t756,t757,t758,t759,t76,t760,t761,t762,t763,t764,t765
           type(dd_complex) ::  t766,t767,t768,t769,t77,t770,t771,t772,t773,t774,t775,t776
           type(dd_complex) ::  t777,t778,t779,t78,t780,t781,t782,t783,t784,t785,t786,t787
           type(dd_complex) ::  t788,t789,t79,t790,t791,t792,t793,t794,t795,t796,t797,t798
           type(dd_complex) ::  t799,t8,t80,t800,t801,t802,t803,t804,t805,t806,t807,t808
           type(dd_complex) ::  t809,t81,t810,t811,t812,t813,t814,t815,t816,t817,t818,t819
           type(dd_complex) ::  t82,t820,t821,t822,t823,t824,t825,t826,t827,t828,t829,t83
           type(dd_complex) ::  t830,t831,t832,t833,t834,t835,t836,t837,t838,t839,t84,t840
           type(dd_complex) ::  t841,t842,t843,t844,t845,t846,t847,t848,t849,t85,t850,t851
           type(dd_complex) ::  t852,t853,t854,t855,t856,t857,t858,t859,t86,t860,t861,t862
           type(dd_complex) ::  t863,t864,t865,t866,t867,t868,t869,t87,t870,t871,t872,t873
           type(dd_complex) ::  t874,t875,t876,t877,t878,t879,t88,t880,t881,t882,t883,t884
           type(dd_complex) ::  t885,t886,t887,t888,t889,t89,t890,t891,t892,t893,t894,t895
           type(dd_complex) ::  t896,t897,t898,t899,t9,t90,t900,t901,t902,t903,t904,t905
           type(dd_complex) ::  t906,t907,t908,t909,t91,t910,t911,t912,t913,t914,t915,t916
           type(dd_complex) ::  t917,t918,t919,t92,t920,t921,t922,t923,t924,t925,t926,t927
           type(dd_complex) ::  t928,t929,t93,t930,t931,t932,t933,t934,t935,t936,t937,t938
           type(dd_complex) ::  t939,t94,t940,t941,t942,t943,t944,t945,t946,t947,t948,t949
           type(dd_complex) ::  t95,t950,t951,t952,t953,t954,t955,t956,t957,t958,t959,t96
           type(dd_complex) ::  t960,t961,t962,t963,t964,t965,t966,t967,t968,t969,t97,t970
           type(dd_complex) ::  t971,t972,t973,t974,t975,t976,t977,t978,t979,t98,t980,t981
           type(dd_complex) ::  t982,t983,t984,t985,t986,t987,t988,t989,t99,t990,t991,t992
           type(dd_complex) ::  t993,t994,t995,t996,t997,t998,t999

           type(dd_complex) :: result

      t1 = ddreal(1) - epinv
      t2 = ddreal(8)
      t3 = (gb)**2
      t4 = (t3)**2
      t5 = (gw)**2
      t6 = ddreal(9) * (t5)**2
      t7 = propZ25 * s25
      t8 = t7 * (t4 - t6) - t2 * t3 * t5
      t9 = intHs16s25s26s34s56x1211D4eps0()
      t10 = -t8
      t11 = s16 + s26 - s34 - s45
      t12 = -s13 + s25 + s26 + s45 + s46 + s56
      t13 = s12 + s16 + s26 - s34 - s45
      t14 = intHLs16s25s26s34s56x1112D4eps0()
      t15 = ddreal(6)
      t16 = ddreal(4)
      t17 = t7 * t15
      t18 = t17 + t16
      t19 = t7 * t3
      t20 = t6 * t7
      t21 = t3 * (t18 * t5 + t19) + t20
      t22 = intHs16s25s26s34s56x1111D2eps0()
      t23 = s12 + s16 + s26 + s46 + s56
      t24 = s12 + s16 + s26 - s34 - s45 + s56
      t25 = s16 + s26 - s34 + s56
      t26 = ddreal(2)
      t27 = -t26 + epinv
      t28 = intHs160s26s34s56x1012D2eps0()
      t29 = intHs160s26s34s56x1021D2eps0()
      t30 = intHs160000x0211D2eps0()
      t31 = intHs160000x0211D2eps1()
      t32 = intHs16s25s26s34s56x1211D4eps1()
      t33 = intHLs16s25s26s34s56x1112D4eps1()
      t34 = t26 + epinv
      t35 = intHLs16s25s26s34s56x1211D2eps1()
      t36 = t1 * s12
      t37 = intHs160s26s34s56x1012D2eps1()
      t38 = intHs16s25s26s34s56x1112D2eps1()
      t39 = epinv * t13
      t40 = s12 - t39
      t41 = t26 * t3 * t5
      t18 = t7 * (t4 + t6) + t41 * t18 / ddreal(2)
      t42 = -s13 + s25 + s26 + s45 + s46
      t43 = t26 * s56
      t44 = -t42 - t43
      t45 = s34 + s45 + s46 + s56
      t46 = t12 * t45
      t47 = s12 * t44 + t46
      t48 = s25 + s26 + s56
      t49 = intHL0s25s26s34s56x1310D4eps1()
      t50 = intHs16s25s26s34s56x1310D4eps0()
      t51 = intHL0s25s260s56x1031D4eps0()
      t52 = intHs16s25s26s34s56x1130D4eps0()
      t53 = intHLs16s25s26s34s56x1111D2eps0()
      t54 = t26 * s45
      t55 = s25 - s12 - s16 + s34 + s56 - s13 + s46 + t54
      t56 = intHs16s25s26s34s56x1141D6eps0()
      t57 = s25 - s12 - s13
      t58 = -s13 + s25
      t59 = s16 - s34
      t60 = t26 * s26
      t61 = ddreal(3)
      t62 = t61 * s46
      t63 = t26 * (s26 + s45 - s13)
      t64 = s16 * t58
      t65 = t58 * s34
      t66 = (s56)**2
      t67 = s56 * t66
      t68 = t26 * t66
      t69 = intHs16s25s26s34s56x1311D4eps1()
      t70 = (s12)**2
      t71 = s12 * t70
      t72 = epinv * t12
      t73 = intHLs16s25s26s34s56x1113D4eps1()
      t74 = s34 + s45
      t75 = t26 * t74
      t76 = s16 + s26 - s56 - s46 - t75
      t77 = t39 * s12
      t78 = intHLs16s25s26s34s56x1112D2eps1()
      t79 = intHs16s25s26s34s56x1211D2eps1()
      t80 = t72 * t45
      t81 = intHLs16s25s26s34s56x1121D2eps0()
      t82 = s16 + s26 + s46
      t83 = s45 + s46 + s56
      t84 = s45 + s46
      t85 = s46 * t84
      t86 = s26 * s34
      t87 = s25 * t45
      t88 = s12 * t48
      t89 = s26 + s45 + s56
      t90 = t84 * s26
      t91 = t89 * s56
      t92 = t91 + t90
      t93 = t62 * s56
      t94 = intHLs16s25s26s34s56x1112D2eps0()
      t95 = intHLs16s25s26s34s56x1113D4eps0()
      t96 = intHs160s26s34s56x1020D2eps0()
      t97 = intHs16s25s26s34s56x1231D6eps0()
      t98 = s12 + s13
      t99 = t16 * s16
      t100 = s12 + s16
      t101 = t61 * s25
      t102 = t26 * s16
      t103 = t2 * s25
      t104 = t61 * s13
      t105 = t16 * s25
      t106 = t61 * s12
      t107 = s16 + s26
      t108 = s26 + s13
      t109 = t61 * s16
      t110 = t15 * s25
      t111 = t26 * s25
      t112 = s13 - t111
      t113 = t26 * s12
      t114 = t2 * s16
      t115 = -s16 + s25
      t116 = ddreal(5)
      t117 = t26 * s13
      t118 = t116 * s26
      t119 = t98 * s13
      t120 = t98 * s16
      t121 = (s26)**2
      t122 = s26 * t121
      t123 = (s16)**2
      t124 = s16 * t123
      t125 = t26 * t123
      t126 = ddreal(7) * t121
      t127 = s34 - s45
      t128 = s25 + s45
      t129 = -s25 + s34
      t130 = s26 * t129
      t131 = s16 + s45
      t132 = s16 + s25
      t133 = s34 * s45
      t134 = (s34)**2
      t135 = s34 * t134
      t136 = t15 * t121
      t137 = s13 + s34
      t138 = ddreal(7) * s26
      t139 = t26 * t128
      t140 = s25 - s16 - s13
      t141 = s26 - s46
      t142 = t61 * s45
      t143 = s26 + s46
      t144 = t61 * s34
      t145 = t15 * s16
      t146 = s16 + s46
      t147 = s45 * s46
      t148 = (s46)**2
      t149 = s46 * t148
      t150 = s16 * s26
      t151 = s34 * s46
      t152 = s16 * s45
      t153 = (s45)**2
      t154 = ddreal(18) * s46
      t155 = ddreal(10) * s45
      t156 = ddreal(11) * s16
      t157 = ddreal(11) * s26
      t158 = ddreal(10) * s46
      t159 = ddreal(7) * s16
      t160 = t116 * s13
      t161 = t16 * s34
      t162 = s12 * s13
      t163 = -(-s13 + s25 - s26 + s34) * s45 + (s13 + t113 + t109) * s25
     & - (-s13 - t111 + t114 + t106) * s26 - (t115 * t61 - s12 + s34 - t
     &117 - t118) * s34 - t119 - t120 - t125 - t126
      t164 = (s13 - t105 + t106) * s26 + (-t103 + t104) * s16 - (s13 - t
     &102 + t101 - s12) * s12 - s13 * t58
      t163 = t16 * t67 + t26 * (-s13 * t123 + t149) - (t120 + t162) * s1
     &3 + (t100 * (s13 + t99) + t70) * s25 - t164 * s26 - (t112 * s34 - 
     &s12 * (t108 - t101) - (s13 - s25 + t109 + t60) * s13 + t110 * t107
     &) * s34 - t163 * s45 - (-(-t137 * t61 + s12 + t138 + t139 + t99) *
     & s46 - t16 * s16 * t128 + t61 * s34 * t132 - ddreal(7) * s26 * t13
     &1 + t116 * t130 + t26 * (-t123 + t133) - s12 * (s16 + s25 - s34 + 
     &s45) - (t127 * t26 + s13 - s25 - t109 - t118) * s13 - t134 - t136)
     & * s46 - (-(s12 + t111 - t161 - t160 + t142 + t159 + t158 + t157) 
     &* s56 - t16 * s25 * t146 - (t155 + t154) * s26 + t116 * (-t147 + t
     &130) - t2 * (t148 + t150) + t61 * (s25 + s16 + s45) * s34 - t125 +
     & ddreal(7) * t151 - ddreal(7) * t152 - ddreal(7) * t121 - (-t141 *
     & t26 - s12 + t140 + t142) * s12 - (-t143 * t2 + s13 - t128 + t144 
     &- t145) * s13 - t134 + t153 - t156 * s46) * s56
      t164 = t26 * t107
      t165 = s12 - s34 - s45 + s56 + s46 + t164
      t166 = intHs160s26s34s56x1021D2eps1()
      t167 = s12 + t102
      t168 = s12 + t164
      t169 = s16 + s26 - s34
      t170 = t169 * s45
      t171 = t13 * s56
      t172 = s16 * t100
      t173 = (s26 + t167) * s26
      t174 = intHs16s25s26s34s56x1120D2eps1()
      t175 = s25 + s26 + s45 + s46 + s56
      t176 = (s13 - t175) * t25
      t177 = epinv * s25
      t178 = intHs160000x0121D2eps0()
      t179 = t11 * (s12 + t12) + t70
      t180 = intHs160000x0112D2eps0()
      t181 = (s12 + t11) * s12
      t182 = t181 - t46
      t183 = intHLs160000x0111D2eps0()
      t184 = t26 + t7
      t185 = t7 * t61 * t5 + t184 * t3
      t186 = -s25 + s16 - s34 - s56 + s13 - s46 - t54
      t187 = (t186 + s12) * s12 + t46
      t188 = intHL0s25s260s56x1012D2eps0()
      t189 = -s12 * t44 - t46
      t190 = intHL0s25s260s56x1013D4eps0()
      t191 = intHL0s25s260s56x1021D2eps0()
      t192 = intHs16s25s26s34s56x1121D4eps0()
      t193 = -s13 - s16 + s25 + s34 + s46
      t194 = s45 + s56
      t195 = t26 * t194
      t196 = t11 * t12
      t197 = intHLs16s25s26s34s56x1121D4eps0()
      t198 = s16 + s26 + s34 + s45
      t199 = s46 + s56
      t200 = t26 * t199
      t201 = t198 + t200
      t202 = s26 + s34 + s45
      t203 = -s12 + s34
      t204 = s16 * t175
      t205 = t203 * s56
      t206 = t202 * s25
      t207 = s25 * s56
      t208 = s26 * s45
      t209 = s25 * s46
      t210 = s26 * t199
      t211 = s45 * t199
      t212 = t16 * s46
      t213 = t212 * s56
      t214 = intHs16s25s26s34s56x1112D4eps0()
      t186 = (t186 + s12) * s12 - t12 * t201
      t215 = intHLs16s25s26s34s56x1211D4eps0()
      t216 = intHs16s25s26s34s56x1121D4eps1()
      t217 = t194 * epinv
      t218 = -s25 - s34 - s46
      t219 = t27 * s13
      t220 = t1 * s16
      t221 = t16 * s56
      t222 = intHLs16s25s26s34s56x1121D4eps1()
      t223 = epinv * t199
      t224 = epinv * t74
      t225 = t27 * t107
      t226 = t26 * (s34 + s45 + t223) + t225 + t224
      t227 = -s25 + s34 - s26
      t228 = epinv * s26
      t229 = epinv * t148
      t230 = epinv * t66
      t231 = s25 * s34
      t232 = t227 * s26
      t233 = s46 + s56
      t234 = s26 * t233
      t235 = t27 * s16
      t236 = intHLs16s25s26s34s56x1211D4eps1()
      t237 = epinv * s45
      t238 = t16 * s26
      t239 = t34 * s12
      t240 = intHs16s25s26s34s56x1112D4eps1()
      t241 = t1 * s34
      t242 = t1 * s45
      t243 = t61 * t199
      t244 = intHLs16s25s26s34s56x1121D2eps1()
      t245 = -s13 + s25 + s26
      t246 = s25 + s26
      t247 = s16 + s25 + s26
      t248 = s25 + s34
      t249 = s16 * t246
      t250 = t111 * s26
      t251 = s12 - s13 + s16
      t252 = -s13 + s26
      t253 = t100 * s13
      t254 = (t251 + t111 + s26) * s26
      t255 = s12 + s16 + s25 + s26 - s34
      t256 = t61 * t246
      t257 = s25 + s26 + s45 + s46
      t258 = s13 * s16
      t259 = t246 * s34
      t260 = s12 * s25
      t261 = t247 * s45
      t262 = t255 * s46
      t263 = t245 * s25
      t264 = t263 * s34
      t265 = (-s13 * t247 + s25 * t248 + t121 + t249 + t250 + t261) * s4
     &5
      t266 = ((s25 + t251) * s25 - s34 * t252 + (t247 * t26 - t203) * s4
     &5 - t253 + t254 + t262) * s46
      t267 = ((t257 * t26 + s56 + t251) * s56 + (s34 + t251 + s25) * s25
     & + (s45 - s13 + t102 + t256) * s45 + (t26 * (s12 + s16 + s45) - t1
     &37 + t256 + s46) * s46 - t253 + t254) * s56
      t268 = epinv * s12
      t269 = intHLs16s25s26s34s56x1141D6eps1()
      t270 = s13 - t246
      t271 = t128 * s34
      t272 = s12 * (-s13 + s25 + s26 + s45)
      t273 = s12 + s16 - s34
      t274 = (-t58 - t60) * s34
      t275 = t26 * t273
      t276 = t245 * s34
      t277 = t276 * s26
      t249 = (-s16 * t270 + t152 - t86) * s45 + (-s13 * t59 - t26 * (-t1
     &52 + t86) + t249 - t271 + t272) * s46 + ((-s34 + t102) * s45 + t10
     &0 * t246 - t253 + t274 + t275 * s46) * s56 + t148 * t273 + t66 * t
     &273 - t277
      t278 = intHLs16s25s26s34s56x1132D6eps1()
      t279 = -s13 + s16 + s25
      t280 = t107 + t101
      t281 = s26 + s34
      t282 = s26 * t107
      t283 = s12 - s13
      t284 = t26 * t283
      t285 = t26 * t252
      t286 = s12 + s26 - s34
      t287 = (s25)**2
      t288 = -t255 * t26 + s45
      t289 = t61 * s26
      t290 = s45 - s13
      t291 = t116 * s25
      t292 = t61 * t100
      t293 = t26 * s34
      t294 = t101 * s34
      t295 = s26 + s45 + s46
      t296 = t26 * t295
      t297 = s16 - s26
      t298 = -t26 * t203
      t299 = s25 * t100
      t300 = s26 - s45
      t301 = s25 + s12 - s13
      t302 = t16 * t59
      t303 = t26 * t59
      t304 = (t279 + s26) * s26
      t305 = t304 + t64
      t274 = s26 * t305 + (s45 * t297 - s13 * t297 + s16 * (s25 + t60) -
     & s26 * (s25 + t293)) * s45 + ((t300 + t275) * s46 - t26 * (t253 - 
     &t121 - t274 - t299) + (t58 + t113 + t109) * s26 - (s45 + t245 - t2
     &98 - t109) * s45) * s46 + ((s12 + s26 - s45 + t303) * s56 - t26 * 
     &(-t121 - t274 - t299 + t258) + (t301 + t109) * s26 - (s45 + t245 +
     & t293 - t109) * s45 + (t26 * t300 + t106 + t302) * s46 - t162) * s
     &56 - t277 * t26
      t277 = t280 * s45
      t306 = t288 * s46
      t307 = intHLs16s25s26s34s56x1123D6eps1()
      t308 = -s13 + s25 + s34
      t309 = (t13 + t111) * s46
      t310 = t16 * s45
      t311 = s25 + s46
      t312 = t26 * t143
      t313 = t13 * s46
      t314 = s12 * s56
      t272 = ((s25 - s16 + s34 - s13 + s45) * s45 - t64 + t276) * s26 + 
     &(-t272 - t313) * s46 + t11 * (s46 * (-s25 - s45 + s13 - t60) + (s1
     &3 - t128 - t312 - s56) * s56) - t121 * t279 - t122 - t314 * t311
      t315 = (s12 + s16 + s25 + s26 - s34 - s45) * s46
      t316 = s25 * s45
      t317 = t26 * s25 * (s16 * (-s13 + s25) + (t279 + s26) * s26) - t61
     & * s25 * (t153 + t276) - epinv * t272 - ((s25 + s12 - s16 - s13) *
     & s25 + (s26 + t251) * s26 - (t252 - t111) * s34 + (s12 + s16 - s34
     & + s13 + t101 - s45) * s45 - t253 + t315) * s46 - ((t310 + t144) *
     & s25 + (s25 - s26 + t283 - t102) * s25 + t309 + t207) * s56 - t316
     & * (t308 * t61 + s26 - t102)
      t318 = intHLs16s25s26s34s56x1131D4eps0()
      t319 = s16 + s26 + s46 + s56
      t320 = s26 + s56
      t321 = s34 + s46 + s56
      t322 = s25 * t74
      t323 = s16 * t48
      t324 = s56 * t199
      t325 = s45 * t320
      t326 = t26 * t325
      t327 = s34 * t12
      t328 = t327 * t319
      t329 = (s13 * t320 - s26 * t321 - t322 + t323 - t324 - t326 + t88)
     & * s12 - t328
      t330 = intHs16s25s26s34s56x1123D6eps0()
      t331 = s16 + s26 - s45
      t332 = t331 - t293
      t333 = s25 + s45 + s46
      t334 = s46 * t128
      t335 = s25 * s26
      t336 = t246 * s56
      t337 = s34 - s45 - s46
      t338 = t84 * s56
      t339 = s25 + s45 + s46 + s56
      t340 = t339 * s34
      t341 = s12 * t320
      t342 = intHs16s25s26s34s56x1121D2eps0()
      t343 = s26 + s45
      t344 = s34 * t343
      t345 = t26 * s46
      t346 = t46 * t25
      t347 = intHs16s25s26s34s56x1122D4eps0()
      t348 = t300 - t293
      t349 = s45 - s56 + s34
      t350 = s16 * t348
      t351 = s45 * s56
      t352 = t349 * s34
      t353 = -s34 - s45 - s56
      t354 = t26 * (-s26 * (s45 + t293) + t350 - t351 + t352) + s46 * t3
     &53 + t121 + t123 - t66
      t355 = t107 - t161
      t356 = s34 * t311
      t357 = s34 * t127
      t358 = t175 * s26
      t359 = t116 * s56
      t360 = t15 * s26
      t361 = t26 * t70
      t362 = intHs16s25s26s34s56x1131D4eps0()
      t363 = s16 + s26 + s56
      t364 = t363 - t293
      t365 = s34 * t320
      t366 = s26 * s56
      t367 = t333 * s26
      t368 = s12 * s34
      t369 = t11 * t25
      t370 = intHLs16s25s26s34s56x1122D4eps0()
      t371 = t320 + t111
      t372 = s26 + t111
      t373 = s34 + s46
      t374 = s16 * t372
      t375 = s26 + s45 + s46 + s56
      t376 = t375 * s56
      t377 = s34 + s45 - s56 - s46
      t378 = t377 - t60
      t379 = s26 - s34 - s45
      t380 = s26 + s46 + s56
      t381 = t379 * t380
      t382 = s12 * (s26 + t43)
      t383 = t371 * s12
      t384 = t11 * t319
      t385 = t384 * t175
      t386 = intHs16s25s26s34s56x1132D6eps0()
      t387 = s34 - s56
      t388 = t26 * t82
      t389 = s25 * t320
      t390 = t333 * s46
      t391 = s56 * t320
      t392 = s26 * s46
      t393 = t116 * s46
      t394 = t384 * t12
      t395 = intHs16s25s26s34s56x1121D2eps1()
      t396 = epinv + ddreal(1)
      t397 = epinv * s56
      t398 = t396 * s46
      t399 = t396 * s25
      t400 = t396 * s45
      t401 = t396 * t148
      t402 = t228 * s34
      t403 = t241 * s45
      t404 = intHs16s25s26s34s56x1122D4eps1()
      t405 = epinv * s34
      t406 = t1 * s26
      t407 = t237 + t406
      t408 = t16 * epinv
      t409 = t408 * s34
      t410 = t1 * t123
      t411 = -t26 * epinv * ((-s34 - s45 + s56) * s34 + t351) + s16 * (-
     &t26 * t407 + t377 - t409) - ((t54 + t161) * epinv - t377 + t406) *
     & s26 - (t397 - t241 - t242) * s46 - (-s34 - s45 + t397) * s56 - t4
     &10
      t412 = t61 * s56
      t413 = t27 * s25
      t414 = t26 * t84
      t415 = epinv * t134
      t416 = t380 * s25
      t417 = t1 * s46
      t418 = t1 * s25
      t419 = s46 + t405
      t420 = t402 * t15
      t421 = s12 * (-t116 * s56 * t419 - t16 * ((-s45 - s56 - t418 - t41
     &7) * s34 + t351 + t90) - t26 * ((-s26 + t237) * s34 + t85 - t415 +
     & t416) - t61 * t391 + epinv * (t358 + t351) + s13 * (-t16 * t241 -
     & t225 + t345 + t412) + s16 * (epinv * (t375 - t293) - s56 + t413 -
     & t414) - t420) + t12 * t411 + t361 * (s26 + s56 - t405)
      t422 = intHs16s25s26s34s56x1131D4eps1()
      t423 = t26 * t320
      t424 = t397 + t225
      t425 = t26 * t83
      t426 = -s25 + t397
      t427 = t15 * s46
      t428 = t427 * s56
      t429 = intHLs16s25s26s34s56x1122D4eps1()
      t430 = s25 - s45 - s56
      t431 = -t430 + t345
      t432 = t16 * t74
      t433 = -s26 + t243 + t432
      t434 = s12 + s45 + s56
      t435 = -s45 + s56
      t436 = s26 * t194
      t437 = s25 * t433
      t438 = t435 * s46
      t439 = -t378
      t440 = s25 + s45 + s26
      t441 = s25 + s56
      t442 = t441 * s26
      t443 = t59 * s46
      t444 = t121 * s46
      t445 = t84 * t66
      t446 = t100 * s56
      t447 = s25 * t148
      t448 = t128 * s56
      t449 = s16 * s25
      t450 = t2 * s45
      t451 = ddreal(7) * s46
      t446 = s25 * ((-t451 - t450) * s56 - t427 * s45) - t16 * (s25 * (s
     &45 * t440 + t66) + t231 * (s26 + s45 + s56 + s25) + t442 * s46) + 
     &t26 * (((-s12 - s34) * s46 - t314) * s25 + ((-s12 - s16 + s34 - s4
     &5 - s46) * s46 - t66) * s26 + s46 * (-t443 + t133) + (t151 - t152)
     & * s56 - t444 - t445 - t85 * s12) - t61 * ((t446 + t152 + t287) * 
     &s46 + s56 * t287 + t447 + t448 * s26) - epinv * ((t26 * (t208 + t3
     &76 + t322) + (-s25 + s34 + s46) * s26 - t374 - t383) * s12 + s13 *
     & ((s16 + t439) * s16 + t381 - t382) - t385) + s13 * (t26 * t286 * 
     &s46 + s16 * t431 + s56 * t434 + t436 + t437 + t438) - s16 * (t153 
     &- t287) - (s45 * t131 - t287 + t446 - t449) * s26 - (s45 * t435 + 
     &t449) * s46 - ((s12 + s16 + s56) * s56 + t153) * s56 + t121 * t430
     & - t148 * t435
      t452 = intHLs16s25s26s34s56x1131D4eps1()
      t453 = s16 * t83
      t454 = t194 * s34
      t455 = t147 + t66
      t456 = s34 - s46
      t457 = t343 * s46
      t458 = s46 * t456
      t459 = -s12 - s16 - s26
      t460 = intHs16s25s26s34s56x1310D4eps1()
      t461 = intHL0s25s260s56x1031D4eps1()
      t462 = intHs16s25s26s34s56x1112D2eps0()
      t463 = intHs16s25s26s34s56x1130D4eps1()
      t464 = intHLs16s25s26s34s56x1211D2eps0()
      t465 = intHs16s25s26s34s56x1311D4eps0()
      t466 = intHLs16s25s26s34s56x1114D6eps1()
      t467 = intHs16s25s26s34s56x1211D2eps0()
      t468 = intHs160s26s34s56x1022D4eps0()
      t469 = s12 + s16 + s26 - s45
      t470 = (t107 - t293) * s45
      t471 = (t469 - t293) * s56
      t472 = s34 * t203
      t473 = t61 * t107
      t474 = -t473 * s34 + t26 * t472 + t172 + t173 - t470 + t471
      t475 = intHs16s25s26s34s56x1220D4eps0()
      t476 = -s13 + s16
      t477 = s25 - s16 - s26 + s34
      t478 = t26 * t246
      t479 = t167 * s25
      t480 = (-t252 - t111) * s34
      t481 = t169 * s46
      t482 = (-t476 - t111 - s26) * s26 + s45 * t477 + (-s16 + s34 - s45
     & + s13 - s46 - t478 - s56) * s56 + t258 - t479 - t480 - t481
      t483 = intHLs16s25s26s34s56x1141D6eps0()
      t484 = intHLs16s25s26s34s56x1132D6eps0()
      t485 = intHLs16s25s26s34s56x1123D6eps0()
      t486 = intHs16s25s26s34s56x1321D6eps0()
      t487 = t98 + t102
      t488 = s13 + s25
      t489 = s13 * (s13 - t132)
      t490 = t26 * t100
      t491 = t26 * t476
      t492 = t61 * t121
      t493 = s12 - s16
      t494 = t26 * t58
      t495 = -t129
      t496 = t476 * s25
      t497 = t26 * t121
      t498 = s25 + s16 - s34
      t499 = t26 * t132
      t500 = t61 * t476
      t501 = t2 * s26
      t502 = t26 * t495
      t503 = -s12 + s16 - s13
      t504 = t16 * t148
      t505 = t116 * t143
      t506 = (s13)**2
      t507 = (-s13 + s16 + s25 + s26) * s45
      t508 = intHs16s25s26s34s56x1221D4eps0()
      t509 = s25 - s16 - s26 + s34 - s56
      t510 = -t377 + t164
      t511 = s25 - s16 + s34 + s45
      t512 = t26 * t380
      t513 = s46 * s56
      t514 = s25 * t199
      t515 = s25 + s45 + s56
      t516 = t74 * s46
      t517 = t248 * s26
      t518 = t515 * s34
      t519 = s16 + s26 - s34 - s45 - s46 - s56
      t520 = s56 - s46
      t521 = t116 * s34
      t522 = t26 * t379
      t523 = t300 * s34
      t524 = s34 * s56
      t525 = t212 * s45
      t526 = t15 * s56
      t527 = t84 * s45
      t528 = s34 * t83
      t529 = t16 * t121
      t530 = t116 * t528 - t26 * t324 + (s56 + t414 + t521 - t289 + s25)
     & * s25 + s26 * (-t16 * t84 + t521 - t526) + t527 - t529 - t142 * s
     &56
      t531 = s34 - s26
      t532 = s56 * t531 - t208
      t533 = -s34 - s56
      t534 = t533 * s45
      t535 = t123 * t175
      t536 = s45 * t66
      t537 = -s25 - s26 - s46
      t538 = t335 * s46
      t539 = t134 * s56
      t540 = (-t134 - t148) * s45
      t541 = -t351 + t86
      t542 = t316 * s46
      t543 = t121 * s34
      t544 = t542 - t543
      t545 = t520 * s45
      t546 = -t545 + t209 + t287 - t148
      t547 = t66 + t148
      t548 = s45 - s56 - s25
      t549 = t519 * s13
      t550 = s56 * t547
      t551 = intHs16s25s26s34s56x1312D6eps0()
      t552 = t515 + t60
      t553 = s25 * t379
      t554 = t101 * s16
      t555 = s25 * (s26 - s34 - s45 - s46 - s56)
      t556 = t151 - t351
      t557 = t393 * s26
      t558 = (t427 + t238) * s56
      t559 = t209 * s56
      t560 = t311 * s26
      t561 = t128 * s26
      t562 = t561 * s46
      t563 = s25 - s46
      t564 = t563 * s46
      t565 = t564 + t207
      t566 = s25 * t134
      t567 = s46 * t387
      t568 = s16 * s56
      t569 = s34 * t148
      t570 = s56 * t143
      t571 = t15 * ((t148 + t570) * s56 + t444)
      t572 = ddreal(13) * s46
      t573 = ((t260 + t26 * t553 + (t552 - s13) * s13 + s46 * t515 - t12
     &1 + t148 + t207 - t366 + t554) * s12 + t116 * ((t335 + t457) * s16
     & + (s26 * t440 - t151) * s56 + t562) - t16 * ((-s56 * t440 - t148 
     &- t209) * s16 - t559 + t560 * s34) + t26 * (s25 * t123 + (-s25 * t
     &169 - t524 + t66 - t86) * s45 + t66 * (s25 - s34 + s56) + t149 + t
     &447) - t61 * (((s16 + s56 + s46) * s34 - t121) * s25 + (t567 - t28
     &2) * s45 + (t86 - t568) * s56 + t569) - (-t26 * t556 + t61 * t547 
     &+ s16 * (-s25 + t43 + t62) + (s26 + s34 + s45 + s46) * s45 - t524 
     &- t555 + t557 + t558 + t549) * s13 - s45 * t565 + t121 * t59 - t15
     &3 * (-s25 - s16 + s46) + t122 + t566 + ddreal(7) * s46 * (t392 + t
     &568) + t571 + t572 * t366) * s12 + t196 * t380 * t510
      t574 = intHLs16s25s26s34s56x1213D6eps0()
      t575 = intHLs16s25s26s34s56x1221D4eps0()
      t576 = t116 * t246
      t577 = t16 * t320 + s16 + s25 + s34 + s45 + t345
      t578 = s34 + s45 + s46
      t579 = s16 * (t84 - t576 - t221)
      t580 = t578 * s46
      t581 = s26 * t246
      t582 = t141 * s56
      t583 = t15 * t325
      t584 = t26 * (-t582 + t392)
      t585 = t116 * t259
      t586 = t61 * (t581 - t524)
      t587 = t316 * t16
      t588 = ddreal(9) * s34
      t589 = t61 * t74
      t590 = t84 * s34
      t591 = s56 * t533
      t592 = s25 - s45
      t593 = t531 * s34
      t594 = t199 * s46
      t595 = s26 * t320
      t596 = s45 * t143
      t597 = t15 * s45
      t598 = t597 * s56
      t599 = -t116 * t595 + t16 * t596 + t26 * (-s46 * t281 - t123 + t13
     &3) - s16 * (t373 * t61 + t138 + t526 + t592) + t153 - t555 + t593 
     &+ t594 + t598
      t600 = t86 * t246
      t601 = s26 + s56
      t602 = t601 * s45
      t603 = s26 - s34 + s45
      t604 = t123 * (s25 + s26 - s45 - s46)
      t605 = t316 * t300
      t606 = t435 * s56
      t607 = t83 * t134
      t608 = s16 * (s25 * (t199 + t310 + t588 - t60) + t15 * t338 - t16 
     &* t591 - t26 * (t121 + t153 - t148) + t61 * (s26 * (s56 + t345 + t
     &589) + t590))
      t609 = t86 * t199
      t610 = t366 * t343
      t611 = t311 * s56
      t612 = t147 * s56
      t613 = s26 * t66
      t614 = (t611 + t209) * s34
      t615 = ddreal(7) * t609
      t616 = intHLs16s25s26s34s56x1212D4eps0()
      t617 = s25 - s16 + s34 + s45 + s56 + t60
      t618 = ddreal(7) * s25
      t619 = s34 - s45 - s56
      t620 = t533 * s26
      t621 = t619 * s46
      t622 = (-t617 + s13) * s13
      t623 = s25 * t531
      t624 = t116 * t623
      t625 = t26 * t514
      t626 = t61 * t375
      t627 = t16 * t456
      t628 = t61 * t320
      t629 = t295 * s45
      t630 = t16 * t387
      t631 = ddreal(12) * s26
      t632 = ddreal(10) * s34
      t633 = ddreal(7) * s45
      t634 = t15 * t199
      t635 = t2 * s56
      t636 = ddreal(9) * s46
      t637 = s45 - s46 - s56
      t638 = t86 * s45
      t639 = t147 * t387
      t601 = t601 * s34
      t640 = s25 * t153
      t641 = t248 * s45
      t642 = s26 * t134
      t643 = t199 * t134
      t644 = t456 * t153
      t645 = ddreal(9) * t231
      t646 = t111 * t71
      t542 = ((s16 * (t375 + t618) - t121 + t148 - t454 + t620 - t621 - 
     &t622 - t587 + t625 - t624) * s12 + (t572 * s56 - t645) * s26 + t11
     &6 * (t121 * t441 + t208 * (-s25 + s56 + s46) + t231 * t637) + t15 
     &* ((t148 + t570 - t130) * s56 + t392 * t246) + t16 * s46 * (-t151 
     &+ t207) - t26 * (t641 * s56 - t128 * t66 - t149 - t447 + t542 - t6
     &40 - t67) - t61 * (s34 * (t66 - t231) + t121 * t127 + t638 + t639)
     & - ddreal(7) * s46 * (-t392 + t601) + (t26 * ((s26 - s45 + s56) * 
     &s34 - t351) - t61 * (t123 + t66 - t458) - t557 - t558 + s16 * (s25
     & + t54 - t628 + t627) - t134 + t555 - t629 - t549) * s13 + ((t626 
     &+ t618) * s16 + s25 * (-t633 + t634 - t632 + t631) + s26 * (-t127 
     &* t16 + t635 + t636) - t83 * (s45 - t393 + t630) + t529) * s16 + t
     &122 - t540 + t642 + t643 + t644) * s12 + t394 * t510 + t646
      t558 = intHs16s25s26s34s56x1212D4eps0()
      t647 = s26 - s34 + s46 + s56
      t648 = t592 * s56
      t649 = t592 * s46
      t650 = s34 * t199
      t651 = -t26 * (-s16 * t647 + t650) + (-t637 - t144 + s26) * s26 + 
     &t123 + t134 + t271 + t316 + t648 + t649
      t652 = s25 + s34 + s56
      t653 = t652 * s26
      t654 = s34 * t143
      t655 = t387 * s56
      t656 = (t118 + t526) * s46
      t657 = s34 - s56 - s46
      t658 = ddreal(9) * s56
      t659 = t16 * t581
      t660 = s25 + s45
      t661 = t129 * s45
      t662 = t316 * s26
      t663 = (s34 * t440 + t153) * s34
      t664 = t311 * s46
      t665 = s25 * t121
      t666 = t321 * t287
      t667 = ddreal(13) * s56
      t668 = t12 * t510
      t562 = (-(t622 + t121 + t518 + t621 + t653 - t148 - t204) * s12 + 
     &t116 * (t366 * t440 - t543 + t562 - t569) + t15 * ((t664 + t570) *
     & s56 + t444) - t16 * s34 * (s26 * t660 + t66) + t26 * ((t134 + t31
     &6) * s46 + (t134 - t661) * s56 + t149 + t536 + t67 + t662 + t663) 
     &+ t61 * ((-t151 + t66 + t148 - t524 + t133) * s25 + s45 * (-t567 +
     & t121) + t535) - ddreal(9) * t524 * t143 + (t16 * (-t366 + t654) -
     & t26 * (t134 + t351 + t87) - t61 * (-t655 + t123 + t148 + t133) + 
     &s16 * (t116 * t456 + t128 - t221 - t289) + t335 - t629 - t656 - t5
     &49) * s13 + s16 * (s26 * (-t116 * t127 + t158 + t658) - t83 * (t11
     &6 * t387 - t427) - t291 * t657 + t659) + (t287 + t148) * s45 + t15
     &3 * t563 + t122 + t665 + t666 + t392 * (t451 + t667 - t632)) * s12
     & + t668 * t651
      t669 = intHs16s25s26s34s56x1222D6eps0()
      t670 = t15 * s34
      t671 = s16 + t512
      t672 = ddreal(9) * s25
      t673 = -t531
      t674 = ddreal(7) * t673
      t675 = ddreal(7) * s56
      t676 = -t15 * t458
      t677 = ddreal(17) * s26
      t678 = ddreal(10) * t199
      t679 = ddreal(13) * s34
      t680 = ddreal(16) * s56
      t681 = ddreal(13) * s45
      t682 = ddreal(11) * s45
      t683 = ddreal(9) * t147
      t684 = t116 * t121
      t685 = t231 * s26
      t686 = t148 * s56
      t687 = s26 * t153
      t688 = t153 * s56
      t689 = s34 * t153
      t690 = t153 * s46
      t691 = t506 * t519
      t692 = s45 * t148
      t693 = s26 * t148
      t694 = s56 * t456
      t695 = s25 * t66
      t696 = ddreal(15)
      t697 = ddreal(10) * s26
      t698 = -t588 + t697
      t660 = s26 * (s45 * t698 + ddreal(30) * t513) - t116 * (-t194 * t6
     &6 + t662) + t15 * t133 * (s25 - s46 - s56) + t16 * (t149 + t566 + 
     &t447 + t695 - t543) - t2 * ((t151 - t121 + t694) * s25 + t569) - t
     &26 * (-t316 * t637 - t122 - t607 + t687 + t688 - t689 + t690 + t69
     &1 - t692) - ddreal(14) * s46 * ((s34 - s26) * s26 - t66) + t696 * 
     &((t366 - t151 + t208) * s56 + t693) - ddreal(13) * t685 + ddreal(1
     &3) * t686 + ddreal(11) * s26 * (s46 * t660 + t207) - ddreal(12) * 
     &t366 * t531 + ddreal(7) * s56 * (t147 - t524) - s13 * ((t572 + t14
     &2 + t675) * s56 - t116 * (t365 - t123) + t26 * (t121 + t134 - t555
     & + t527) + s16 * (t199 * t2 - t111 - t142 + t674) + t208 + ddreal(
     &11) * t210 + t676 + t161 * s45) + ((t116 * t83 + t289 + t672) * s1
     &6 + (ddreal(20) * s46 + t682) * s56 - ddreal(7) * t528 + ddreal(10
     &) * t547 + s25 * (-t633 - t679 + t678 + t677) + s26 * (-t161 + t68
     &1 + t154 + t680) - t153 + t683 + t684) * s16 + t642
      t699 = (-t320 + t111) * s12
      t611 = (((t118 - t670) * s25 - t16 * (t595 + t316) + t26 * (t611 +
     & t390 - t506 - t590) + t325 * t61 + s13 * (t26 * t511 + s56 + t289
     &) + s16 * (-s26 + t414 + t103) - t524 + t699) * s12 + t660) * s12 
     &+ t196 * t510 * t671
      t660 = intHs16s25s26s34s56x1212D4eps1()
      t617 = -t617
      t700 = (-ddreal(7) + epinv) * s25
      t701 = t1 * t637
      t702 = t405 * t61
      t703 = t657 * epinv
      t704 = -t123 - t121
      t705 = t704 * t1
      t706 = t405 * t199
      t377 = -t26 * t706 + (t649 + t648 + t134 + t271 + t316) * epinv + 
     &s16 * (-t26 * (t703 + t406) + t377) + s26 * (s34 + t701 - t702) + 
     &t211 + t650 + t705
      t707 = s56 + t406
      t708 = epinv * t456
      t709 = -t237 + t418
      t710 = t26 * epinv
      t711 = -ddreal(1) + t710
      t712 = t1 * t153
      t713 = t392 * t1
      t714 = t405 * t143
      t715 = t366 * t1
      t716 = t351 * t1
      t717 = t513 * t1
      t718 = t405 * s45
      t719 = t405 * s56
      t720 = t1 * t148
      t721 = t1 * t66
      t722 = -t116 * t713
      t723 = -t15 * t717
      t724 = -t61 + epinv
      t725 = t116 * epinv
      t726 = t16 - t725
      t727 = t726 * s34
      t728 = s46 - t397
      t729 = t427 * epinv
      t730 = ddreal(10) * epinv
      t731 = t529 * t1
      t732 = t61 * epinv
      t733 = -ddreal(7) + t732
      t734 = t626 * t1
      t735 = -s34 - t228
      t736 = t1 * t121
      t737 = t208 * t1
      t738 = t735 * s25
      t739 = t316 * s34
      t740 = t242 - t177
      t741 = s34 - t242
      t742 = t177 * t148
      t743 = t406 + t417
      t744 = s25 + t406
      t745 = t177 * s46
      t746 = t743 * s56
      t747 = -ddreal(1) + epinv
      t748 = t316 * t396
      t749 = t1 * t149
      t750 = t1 * t67
      t751 = t1 * s13
      t752 = t751 * t519
      t753 = t1 * t122
      t666 = (-t16 * (-t715 - t714) - t26 * ((s26 - s45 + s56 + t405) * 
     &s34 - t716) - t61 * (t151 - t410 - t720 - t719 + t718 - t721) + s1
     &6 * (t116 * t708 - t16 * (s34 - s46 + t397) + t61 * t707 - t54 - t
     &709) - (t45 * t711 + t406) * s25 - (-t406 - t417) * s45 + t134 + t
     &712 - t722 - t723 + t752) * s13 + ((s25 * t733 - t734) * s16 + s25
     & * (t116 * (-s34 * t27 + t223) + t633 - t634 + t238 * t724) + s26 
     &* (-t726 * s45 + t730 * s46 - t635 + t727 - ddreal(9) * t728) + t8
     &3 * (-t116 * t728 + s45 - t221 + t727 + t729) - t731) * s16 + (epi
     &nv * t287 - t134 - t720) * s45 + t153 * (-s34 + t417 + t177) - t64
     &2 - t643 - t753 + t666 * epinv + t177 * t121
      t726 = ddreal(13) * t715 * s46
      t130 = (-(-(-t121 - t518 - t621 - t653 + t148) * epinv - t1 * ((t6
     &17 + s13) * s13 - s16 * t375) - t121 + t148 - t454 + t620 - t621 -
     & t700 * s16 - t587 + t625 - t624) * s12 - t116 * ((s26 * (s25 + t4
     &05) - t316) * s26 + (t738 + t737) * s46 + (t738 + t737 + t736) * s
     &56 + t739 + t405 * t148) - t15 * ((t746 + t720 - t130 - t745) * s5
     &6 + t392 * t744) - t16 * (t405 * (t66 + t561) - t569 + t559) + t26
     & * ((t663 + t662) * epinv + (t748 + t415) * s46 + (t748 + t403 + t
     &415) * s56 + t66 * (-s25 - t242) - t447 - t640 - t749 - t750) + t6
     &1 * ((-t177 * t435 - t231 + t66) * s34 + (s34 * t740 - t716) * s46
     & + t121 * t741 + t638 + t742 + t177 * t66) + ddreal(9) * s34 * (t3
     &97 * (-s26 - s46) + t335) + ddreal(7) * s46 * (t392 * t747 + t601)
     & + t666 - t726 - t158 * t402) * s12 + t668 * t377 - t646
      t559 = intHLs16s25s26s34s56x1321D6eps0()
      t601 = s25 + s34 + s45 + s46 + t43 + t289
      t625 = (-s25 - s56 + s46) * s26
      t646 = s25 - s45 - s46
      t653 = t646 + t161 + t289
      t663 = t66 + t148 + t133
      t666 = t26 * t663
      t727 = ddreal(7) * s34
      t738 = t524 * t116
      t748 = t320 * t153
      t754 = s46 * t66
      t755 = -t134 - t513
      t756 = t123 * t48
      t757 = s56 * t343
      t758 = t366 * s45
      t759 = t435 * t121
      t760 = s16 * (-t16 * (t621 + t133) + t26 * (t121 + t153 + t66 + t1
     &48) + s25 * (-s45 + t512 - t727) + s26 * (s46 - t54 + t221 - t727)
     & - t738)
      t761 = ddreal(10) * t609
      t762 = t16 * (-s34 * t455 + t392 * t194 + t612)
      t763 = t83 - t293
      t764 = ddreal(9) * t199
      t765 = t26 * (s26 * (t83 - t521) + t134 + t513) + s34 * (s45 - t76
     &4) + t148 - t153 + t66
      t766 = t144 * t380
      t767 = t766 * t439
      t768 = intHs16s25s26s34s56x1213D6eps0()
      t769 = s25 * t300
      t770 = t194 * s46
      t771 = t16 * (-t366 + t151)
      t772 = t116 * s45
      t773 = -t148 + t231
      t774 = t121 + t133
      t775 = t121 * s45
      t776 = s34 * t66
      t777 = -t664 - t207
      t778 = t121 * s56
      t779 = -s25 - s34 + s46
      t780 = t2 * s34
      t781 = t116 * t199
      t782 = t348 * t380
      t439 = t782 * t439
      t783 = t26 * t124
      t784 = ((-s45 + t360 + t781 - t144) * s16 + (-t16 * (t26 * t657 + 
     &s45) + t360) * s26 + s34 * (s34 + s45 - ddreal(7) * t199) - t200 *
     & (-t199 + t54)) * s16 + t439 + t783
      t564 = ((t260 + t26 * t769 - t294 + s16 * (t375 + t105) - s34 * t3
     &75 - t121 + t148 - t366 + t514 - t622 + t770) * s12 - ddreal(7) * 
     &s26 * t773 + t116 * ((t514 + t211) * s26 + s34 * t777 + t778) - t1
     &6 * ((t133 - t207) * s46 + t543 + t638 + t776) + t26 * ((-t335 + t
     &66) * s45 + t441 * t66 + t149 + t447 + t566) + t61 * (s25 * t774 -
     & t351 * t456 + t775) + t571 - ddreal(9) * t151 * t320 - (-t26 * t5
     &34 - t61 * (t655 - t123 - t148 + t86) + t656 - s16 * (t128 - t393 
     &+ t630 - t289) + t134 - t555 + t629 - t771 + t549) * s13 - (-(t291
     & + t626) * s16 - s25 * (t634 - t142 - t727 + t501) - s26 * (-t161 
     &+ t658 + t158 + t772) + t83 * (-t427 - t359 + t161) - t529) * s16 
     &- (-t134 + t564 + t207) * s45 - t153 * t779 + t122 + t642 + t643 +
     & t366 * (t572 - t780)) * s12 + t12 * t784
      t571 = intHLs16s25s26s34s56x1312D6eps0()
      t622 = intHLs16s25s26s34s56x1231D6eps0()
      t626 = t515 + t289
      t629 = t129 * s56
      t630 = t61 * (-t325 + t323)
      t638 = s34 - t425
      t656 = t86 * t2
      t785 = t26 * t638
      t786 = s26 - s56
      t787 = t61 * t786
      t788 = t248 + t787 - t54 - t212
      t789 = t74 * s45
      t790 = t127 * s56
      t791 = ddreal(10) * s56
      t792 = t61 * t153
      t793 = t531 * s45
      t794 = t531 * t148
      t795 = t524 * s45
      t796 = -t134 - t153
      t797 = s16 * ((t791 + t450) * s46 + t116 * t547 + t26 * (-t790 + t
     &121) + s25 * (t781 - t161 + t60) - s26 * (s45 + t627 - t675) - t59
     &0 + t792)
      t798 = t15 * t513 * t199
      t609 = t2 * t609
      t799 = intHLs16s25s26s34s56x1222D6eps0()
      t800 = t139 + t118 + t412
      t801 = s16 * (t43 + t110 + t289)
      t802 = t61 * (t335 - t325)
      t803 = t2 * s46
      t804 = ddreal(11) * s56
      t805 = ddreal(9) * s26
      t806 = ddreal(12) * s46
      t807 = s26 * (-t161 + t142 + t791 + t806)
      t808 = ddreal(16) * s46
      t809 = (t681 + t158) * s46
      t633 = (t808 + t633 + t526) * s56
      t810 = s25 - s56
      t811 = s45 - s46
      t812 = t811 * s46
      t813 = -t513 + t86
      t814 = -t392 + t524
      t815 = ddreal(22)
      t816 = ddreal(11) * s46
      t817 = t513 * t2
      t818 = ddreal(10) * t121
      t439 = t26 * (t124 + t439) + s16 * (s16 * (s45 + t764 + t501 - t14
     &4) + t16 * t547 + s26 * (t199 * t696 - t679 - t772) + t134 - t153 
     &- ddreal(7) * t211 + t817 + t818 - t679 * t199)
      t819 = -t16 * t322 + t26 * (s46 * t339 + t207 - t506 - t86) + s13 
     &* t800 - t524 + t801 + t802 + t383
      t91 = (t819 * s12 + t116 * s56 * ((s26 + s46) * s45 - t524) + t16 
     &* ((t121 + t66 + t148) * s25 + t149 + t133 * t810) - t2 * ((t151 +
     & t271 - t207) * s46 + s56 * (-t121 + t231)) - t26 * ((-s45 * t637 
     &- t134) * s25 + s45 * (t812 - t351) - t122 - t687 + t691 - t775) -
     & t61 * (t208 * t248 - t536 - t67) - ddreal(9) * s26 * ((-s25 - s45
     &) * s46 - t207 + t231 - t66) - ddreal(10) * s56 * t813 + ddreal(10
     &) * t444 - ddreal(13) * s46 * t814 - s13 * ((t805 + t804) * s46 + 
     &t116 * t91 + t26 * (s45 * t578 - t121 - t555) - t61 * (-t208 + t52
     &4) + t676 - (-s16 + s26 + s34 - s45 + t111 - t221 - t803) * s16 - 
     &t86) + ((t83 + t291 + t289) * s16 + s25 * (-t142 + t805 - t727 + t
     &678) - t528 + t807 + t809 + t633 + t792 + t684) * s16 + t642 + t51
     &3 * (s26 * t815 + t816) - t86 * (t360 + t806)) * s12 + t12 * t439
      t676 = s26 - s34 + s46
      t678 = s46 * t531
      t691 = t245 * t107
      t819 = -s16 + s25 + s34
      t820 = s25 + s26 + s34 - s13
      t821 = t820 * t169
      t822 = -s12 - s13 + s25 + s26 + s34
      t823 = t822 + t54
      t824 = s13 + t106
      t825 = -s16 + t113
      t826 = t116 * s12
      t827 = s25 - s34 - s13
      t828 = s12 * t283
      t829 = t824 * s16
      t830 = (t131 * t26 - s12 - s13 + s25 - s34 + s46 + t289) * s46
      t831 = t169 * t153
      t832 = t169 * t148
      t833 = t481 * t823
      t834 = propW16 * t48
      t835 = t834 * t25 * t187 - s12 * t121
      t836 = t113 + t109
      t837 = t61 * t132
      t838 = s16 - s26 - s34
      t839 = t838 - t111
      t840 = t16 * s12
      t841 = ddreal(7) * s12
      t842 = t107 + t840
      t843 = s12 - s45 - s46
      t844 = t20 * (-t26 * t835 - (-s26 * (t283 + t499) + t258 - t479) *
     & s12 - (-s12 * (s13 - t478) + t691 - t276) * s34 - (t819 * s12 + t
     &821) * s45 - ((t279 + t296 - t106 + s56) * s56 - s25 * t825 + (t27
     &9 - t826 + s26) * s26 + (t107 + t106 - s34) * s34 + (t827 + t102 +
     & t289 + s45) * s45 - t828 - t829 + t830) * s56 - t831 - t832 - t83
     &3)
      t845 = intHs16s25s26s34s56x1221D4eps1()
      t846 = t396 * s56
      t847 = epinv * s46
      t848 = t1 * s56
      t849 = (-t848 + t847) * s26
      t850 = t396 * s34
      t851 = t177 * s34
      t852 = epinv * t121
      t853 = s56 * t194
      t854 = -t26 + t732
      t855 = epinv * t128
      t856 = t854 * s16
      t857 = s26 - s45 + s46
      t858 = t228 * s56
      t859 = t857 * s34
      t860 = t300 * epinv
      t861 = t351 * epinv
      t862 = t513 * epinv
      t863 = t392 + t66
      t864 = t435 * s25
      t865 = epinv * t122
      t866 = t852 * s46
      t867 = epinv * t67
      t868 = t134 * s45
      t869 = t134 * s46
      t870 = t1 * t506
      t519 = t870 * t519
      t528 = (-(-t61 * (epinv * t375 + t413) + t425) * s16 - t16 * t85 -
     & t2 * t436 + t528 * t61 - ddreal(7) * t392 - ddreal(7) * t853 - ep
     &inv * t530 - s25 * (-t597 + t781 - t588 + t631) - ddreal(11) * t51
     &3) * s16 - t539 - t689 - s13 * (-t116 * t863 - t2 * t570 + t26 * (
     &-t147 + t859 + t230 + t858 + t415) + t61 * ((s56 - t860) * s34 - t
     &148 + t861 + t862) + (epinv * t520 - t26 * t855 + t61 * (s34 - s46
     & - t406) + s25 + s45 - t526 - t405 * t116 + t856) * s16 - s25 * (s
     &26 * t711 - t224 * t26 + s34 + s45 + t417 + t848) + (-s56 - t406) 
     &* s45 - t121 - t134 + t229 - t712 - t719 + t525 * epinv) + ((-t864
     & - t287) * epinv + t153) * s26 + ((s56 + s46 + t177) * s25 + s45 *
     & (s56 + t177)) * s45 + t865 + t866 + t405 * t546 + t397 * t148 + t
     &867 - t868 - t869 + t519
      t871 = -t151 + t121
      t872 = t231 * s46
      t873 = t397 * t134
      t874 = ddreal(18) * s26
      t699 = (t360 - t521 - t310) * s25 + t26 * ((t846 + t398) * s25 + t
     &849 + t513 * t396) - t61 * (-t325 + t121) + (-t26 * t380 * epinv +
     & s26 - s56 - t220 + t241 + t242 + t418 - t751) * s13 + s16 * (-s26
     & * t34 + t242 + t417 - t700 - t846) + (t850 + t177) * s26 + (t400 
     &- t241) * s46 + t396 * t66 + t401 - t403 + t719 + t851 + t852 + t6
     &99
      t700 = intHLs16s25s26s34s56x1321D6eps1()
      t875 = epinv * t48
      t876 = t405 * (s26 + s46 + s56)
      t877 = t379 * s34
      t878 = s26 - s34 + s56
      t879 = t878 * s46
      t880 = t26 * t877
      t881 = t61 * t876
      t882 = t220 * t763 + t325 + t391 + t87 + t879 - t880 + t881
      t883 = t321 * s56
      t884 = s16 * s34
      t885 = t86 * t16
      t886 = -s56 - t405
      t887 = t147 * t1
      t888 = s26 * t456
      t889 = -t148 - t86
      t890 = -s16 + s45
      t891 = -s26 - s56
      t892 = s16 * s46
      t126 = -epinv * (s34 * ((-t393 - t658) * s46 - t126) - t26 * (t133
     & * t320 - t148 * t343 - t343 * t66 - t642 - t756) + t61 * (t686 + 
     &t748 + t754) - t761 + t762 + s25 * ((-t781 + t142) * s34 - t26 * t
     &755 + s26 * (t199 - t727) + t147 + t148 + t153 + t351 + t66) + (t1
     &53 + t148) * s46 + t607 + t67 + t689 + t758 + t759 + t760) + (-t11
     &6 * t338 + t16 * ((-s26 - s34 + t397) * s45 - t847 * t387) + t26 *
     & (epinv * t663 + t458 + t888) + t61 * (s56 * t886 - t228 * t127 - 
     &t134 + t736 - t887) + (-s56 + t418 + t847) * s26 + t415 - t712 + t
     &413 * t45 + t220 * t653 - t752) * s13 + (t66 - t133) * s16 + (t443
     & + t153 + t568) * s25 - t127 * t134 + t287 * t45 + t149 - t778 - t
     &150 * t435 + t147 * t74
      t126 = (((s34 + t875) * s12 + t26 * (t151 + t121 - t134 + t560) + 
     &t61 * (t366 + t884) + epinv * (-t61 * (t325 - t323 + t86) + s25 * 
     &(s26 - t54 - t144) + s26 * t520 + t121 - t590 - t524 * t26) + s25 
     &* t321 + s45 * t878 + t870 + t883 - t751 * t601 + t885) * s12 - t1
     &16 * ((-t892 - t121) * s34 + t445 + t662 + t207 * t84) + t15 * s26
     & * t755 - t16 * ((s45 * t132 + t134 + t513) * s46 + t231 * t890) -
     & t2 * s26 * (-t884 + t351) + t26 * ((t335 + t121 - t153 - t148 - t
     &316) * s16 + ((s26 + s34 - s46) * s26 - t148) * s25 + s34 * t123 +
     & s45 * t889 + (-t335 - t133) * s56 + t122 - t67 - t689) + t61 * ((
     &s56 * t82 + t148) * s34 + (-t121 - t568) * s45 + (-t892 - t336) * 
     &s56 + t134 * (-s25 - s16 - s56) - t693 - t748) + ddreal(7) * s46 *
     & (s45 * t891 + t86) - t126) * s12 - t668 * t882
      t443 = intHLs16s25s26s34s56x1221D4eps1()
      t445 = s26 + s34 + s56
      t560 = t445 * s45
      t663 = (s26 + s56 - t241) * s46
      t893 = s26 - s45 - s46
      t894 = t893 * s46
      t895 = t257 * s34
      t896 = t26 * t48 * s26
      t897 = t116 * (t66 + t351 + t90)
      t898 = t16 * t199
      t899 = t127 * s46
      t900 = (-t151 + t335 + t448) * s56 + t538
      t901 = t151 - t66
      t902 = s16 + s56
      t903 = s26 * t901 - t147 * t902 - t686
      t904 = t208 * t199
      t905 = t128 * s45
      t257 = (-t392 - t905) * s16 + (t595 - t153) * s34 - t134 * t257 - 
     &t444
      t595 = -s16 * t66 + (s16 * t657 - t148) * s25 + (t884 - t121) * s4
     &5 - t67 - t748 + t884 * t143
      t906 = (s25 + s16 - s56) * t121
      t568 = (-s25 * t83 - t153) * s25 + (-t287 + t334 + t568) * s34 - t
     &148 * t84 + t122 - t539 + t906 + t150 * t430
      t907 = t696 * s26
      t908 = ddreal(9) * s16
      t909 = ddreal(7) * t316 * t320 + ddreal(7) * t693
      t910 = ddreal(10) * t513 * t441 + ddreal(10) * t904
      t911 = t513 * (t907 + t908)
      t257 = t116 * ((t148 + t351) * s16 - t569 + t316 * t373) + t15 * t
     &900 - t16 * t595 - t2 * t903 - t26 * t568 - t61 * t257 + epinv * (
     &(t428 + t656) * s26 - t116 * s45 * (t602 - t121) + t16 * ((t366 + 
     &t899) * s56 + t444) - t26 * ((t812 + t133) * s26 + t66 * t811 + t6
     &04 - t147 * t337) + t61 * (t610 - t642 + t569) + t615 + s25 * ((t8
     &98 - t772 - t144) * s34 + s26 * (t199 + t54 + t780) - s45 * (t199 
     &+ t54)) - s56 * (t134 + t591) + t608 + t686 - t689 - t869 + t133 *
     & t533) - (t16 * (s34 * t811 + t366) - t26 * (-t147 + t86) + t61 * 
     &(s25 * t194 - t121 + t134 + t148 + t209 + t231) + t817 - epinv * t
     &599 - s16 * (s25 + t787 - t54 + t627) + t153 - t335 - t524 + t897 
     &+ t752) * s13 - s34 * (t66 - t351 + t561) + t909 + t910 + t911
      t561 = intHLs16s25s26s34s56x1212D4eps1()
      t568 = t724 * s34
      t591 = t60 * t1
      t595 = -t417 - t848
      t627 = s34 * t74
      t787 = t26 * (t87 + t627) + s16 * (-t701 + t568 + t591) + (s26 * t
     &724 - t417 - t848) * s34 + t595 * s45 + t410 + t736 - t406 * t637
      t812 = t177 * t199
      t617 = t751 * t617
      t817 = ddreal(7) * epinv
      t900 = t61 - t817
      t903 = t406 + t848
      t912 = t417 + t405
      t724 = t724 * t45
      t913 = t241 * s46
      t914 = -t15 + t408
      t915 = t914 * s34
      t916 = -ddreal(1) + t732
      t917 = s34 + t223
      t918 = -s34 - t847
      t919 = t228 * s45
      t920 = t418 + t405
      t921 = t86 * (s46 - t177)
      t922 = ddreal(1) - epinv
      t923 = t747 * s46
      t924 = s34 + s45 - s46
      t644 = (-t540 + t122 + t644 + t643 + t642) * epinv + (-t16 * (-t71
     &5 + t523) - t26 * ((t241 - t242) * s56 - t405 * t300) + t61 * (t13
     &4 - t913 + t410 + t720 + t721) - t722 - t723 - s16 * (-t16 * t912 
     &+ t242 * t26 - t61 * t903 + t418 + t670) - s25 * (t406 + t724) - t
     &415 + t712 + t737 + t887 + t752) * s13 + (-(s25 * t900 + t734) * s
     &16 + s25 * ((ddreal(1) - t817) * s45 + t15 * t917 - t16 * (-s26 * 
     &t916 + s46 + s56) - ddreal(10) * t405) + s26 * (-t914 * s34 - t1 *
     & (t310 + t636 + t635)) - t83 * (-t1 * (s45 - t393 - t221) + t915) 
     &- t731) * s16 + (-t121 + t650) * s25 - t122 + t147 * t924
      t731 = s25 + t237
      t734 = t740 * t66
      t644 = t116 * ((t403 - t852 + t919) * s25 + (t737 + t851) * s46 + 
     &((t406 + t242) * s26 + t851) * s56 - t543) + t15 * ((s26 * t920 + 
     &t720 + t746) * s56 + t392 * (t406 + t418)) + t16 * (((s45 + s56 - 
     &t847) * s56 + t85) * s25 - t241 * t148) - t2 * s56 * (-t209 + t86)
     & + t26 * (s25 * ((s25 + t242) * s45 + t231) + (t641 * epinv + t287
     &) * s56 - t742 + t749 + t750 + t734 + t209 * t731) - t61 * (((s25 
     &- s34 + t242) * s34 - t316) * s26 + s34 * (-t241 * s25 + s45 * t91
     &8 - t151 - t153) + (-t134 - t887) * s56 + t121 * (-t242 - t405) + 
     &t241 * t66) + t726 - ddreal(9) * t921 - ddreal(7) * s46 * ((-t405 
     &+ t923) * s26 + t524 * t922) - t644
      t925 = t177 * t26
      t926 = intHLs16s25s26s34s56x1312D6eps1()
      t927 = t914 / ddreal(2)
      t928 = t129 * s46
      t929 = t233 * epinv
      t701 = t26 * s34 * (t929 + s34 + s45) + (-t26 * (t223 - t406) - s4
     &5 + s46 + s56 + t568 + t220) * s16 + (s34 * t927 - t701) * s26 + (
     &s25 - t417 - t848) * s45 + t231 - t629 + t736 - t928
      t930 = -ddreal(1) + t408
      t931 = -t16 + epinv
      t932 = t83 * epinv
      t933 = t177 * t300
      t934 = t387 * epinv
      t935 = t16 * t711
      t936 = t142 * epinv
      t937 = t725 * s45
      t938 = ddreal(11) * epinv
      t939 = s25 - t237
      t940 = -s34 - t397
      t941 = -s34 - t237
      t942 = t177 * s56
      t943 = t405 * t121
      t944 = t417 + t848
      t945 = t664 + t207
      t946 = t848 * t121
      t947 = t177 * t134
      t948 = t524 * s46
      t949 = s34 - t237
      t740 = t740 * t121
      t746 = ((s25 + t417) * s46 + t746) * s56
      t950 = t779 * epinv
      t951 = t563 * s45
      t733 = -(t61 * t932 + t521 - t898 + (-t61 + t725) * s25 + s26 * t7
     &33) * s16 - t15 * t85 * epinv + ((ddreal(9) - t938) * s46 + t54) *
     & s56 - t116 * (t397 * t194 - t148) + t16 * ((s34 + t932) * s34 + t
     &66) - t2 * t650 - s25 * (t15 * t223 - t898 + s26 * t935 + (t15 - t
     &817) * s34 + s45 - t936) - (t914 * s26 - (-ddreal(13) + t408) * s3
     &4 - t1 * (t658 + t158) - s45 + t937) * s26 - t133 + t62 * s45
      t707 = t733 * s16 + (t134 + t951 + t148) * s45 + t199 * t287 - t51
     &9 + s13 * (t16 * ((s26 - s45 - t847) * s34 - t715) + t26 * (s56 * 
     &t741 + t718) - t61 * (-s34 * t735 - t151 + t410 + t719 + t720 + t7
     &21) + t722 + t723 + s16 * (-t16 * (s46 + t934) - t61 * t707 + t54 
     &+ t670 + t709 + t393 * epinv) + s25 * (t27 * t45 + t406) + t415 - 
     &t712 - t737 - t887) + s34 * (-s34 * t917 + t287) + (t287 + t812 + 
     &t524 - t229 - t415) * s45 + t665 + t753 + t950 * t153 - t228 * t13
     &4
      t709 = t397 * t2
      t707 = ((s16 - s34 - t177) * s12 - t26 * (-s34 * t657 + t933) + t6
     &1 * (t123 + t851) - t885 - s16 * (s25 * t930 + s26 * t931 + s45 - 
     &t200 + t521 + t932) - (s25 + t848 - t405) * s26 - (-epinv * t129 -
     & t242 - t848) * s46 - t231 - t617 + t720 - t736 - t870 + t133 * t3
     &96 + t397 * t129) * s12 + t116 * ((s25 * t944 + s45 * t944) * s26 
     &+ t946 + t405 * t945) + t15 * (s26 * (t134 + t713) + t746) - t16 *
     & (s34 * (s56 * t940 - t316) + (s34 * t941 + t942) * s46 + t569 - t
     &943 + t86 * t939) + t26 * (-t207 * t127 - t209 * t127 + t662 * t39
     &6 + t124 + t689 + t734 - t742 + t749 + t750 - t947) + t61 * (s25 *
     & (s34 * t949 + t148) - t129 * t66 + t740 + t351 * t912) + t726 - d
     &dreal(9) * s34 * (t847 * t891 + t366) - ddreal(7) * (-t720 + t86 -
     & t851) * s26 - ddreal(7) * t948 + t707 + t86 * (t709 - t158)
      t733 = intHLs16s25s26s34s56x1222D6eps1()
      t734 = t127 * epinv
      t737 = t61 * (s34 - s46 + t228)
      t891 = t300 * s46
      t914 = t678 + t86
      t433 = -s25 * t433
      t917 = t657 * s56
      t952 = (-t16 * t223 - s25 + t220 + t423 + t734 - t737) * s16
      t953 = t61 * t914
      t954 = t26 * ((-s26 * t300 - s56 * t300 - t891) * epinv + t133 + t
     &134 - t147 + t366)
      t876 = t16 * t876
      t955 = t366 * t16
      t956 = t27 * s26
      t957 = s26 + t237
      t958 = t957 * s26
      t959 = ddreal(9) * epinv
      t960 = (-ddreal(10) + t959) * s26
      t961 = s45 + t418
      t962 = epinv * t153
      t963 = t242 * t148
      t964 = t402 * s56
      t965 = t335 * t1
      t966 = t311 * s45
      t967 = t406 * t66
      t968 = -s56 + t847
      t633 = s13 * ((-t161 - t675) * s45 + ((-ddreal(12) + t938) * s56 +
     & t960) * s46 + t116 * ((t397 + t237 - t406) * s56 + t151 - t208) -
     & t15 * (s46 * t912 + t66) + t26 * (-t134 - t887 + t718 - t852 - t7
     &12) + t61 * (t958 - t719) - (t927 * s25 - (epinv * t2 - ddreal(7))
     & * s46 - (-t26 + t408) * s56 - t144 + t734 + t956 + t220) * s16 + 
     &s25 * (-t927 * s26 + t26 * (t932 + t568) - t597 - t781) - t402) + 
     &((-(t291 + t289) * epinv + s25 + s26 + t242 + t417 - t846) * s16 +
     & (t816 + t238 + t675) * s45 + t26 * (t366 + t153) - t61 * (s34 * t
     &295 - t66) - epinv * (-s34 * t83 + t633 + t684 + t792 + t807 + t80
     &9) - (s25 + ddreal(10) * t223 + (ddreal(1) + t959) * s26 + t900 * 
     &s34 - s45 - t359 - t427 - t936) * s25 + ddreal(9) * s46 * t143 + d
     &dreal(12) * t513) * s16 + (-t287 - t415) * s26 - t122 - t776 - t87
     &2 + t524 * (s25 - s45 + s34)
      t781 = t248 * epinv
      t807 = t287 * s46
      t809 = t237 * t66
      t900 = t233 * s56
      t936 = t351 * t143
      t969 = t147 * s34
      t970 = epinv * t815
      t971 = (s25 * t696 + (ddreal(23) - t970) * s26 - ddreal(11) * t847
     &) * s56 + ddreal(14) * t392
      t519 = s46 * t971 - t116 * ((t936 - t776) * epinv + t969) + t15 * 
     &((t121 + t231) * s45 + t943) - t16 * (s25 * (-t231 + t852) + (t405
     & * t810 - t287) * s45 + ((-s56 + t177) * s56 - t121) * s56 - t640 
     &+ t742 - t749) - t2 * (((s56 - t847 - t855) * s46 - t942) * s34 - 
     &t662 + t207 * t968 + t397 * t121) + t26 * (((s25 + s34 - s45 - s56
     &) * s34 - t962) * s26 + (s45 * (-t242 + t177) + t134) * s46 + t121
     & * (-s25 - s34 - t237) + t134 * t961 - t519 + t689 - t865 + t963 -
     & t177 * t153 + t861 * t592) + t61 * ((t153 + t287) * s56 - t867 + 
     &t807 - t809 + t208 * (s45 + t781)) + ddreal(10) * ((s25 - t228) * 
     &s26 - t230) * s46 + ddreal(10) * t964 - ddreal(9) * (-t966 - t965)
     & * s56 - ddreal(9) * t921 + ddreal(9) * t967 - ddreal(9) * t392 * 
     &(-s26 + t855) + ddreal(7) * (t316 - t928) * s46 + ddreal(7) * t536
     & + ddreal(13) * t847 * t814 + ddreal(13) * t904 + t633 + ddreal(12
     &) * s46 * (t900 + t402)
      t519 = ((t26 * (-t870 + t85) - t581 * t61 - epinv * (-t26 * (-s25 
     &* (t199 - t75) - s46 * t83 + t86) + t802 + t801 - t524) - (s16 - s
     &34 + s45 - s46) * s56 - t66 - t90 + t751 * t800 - t955 - t268 * t3
     &71) * s12 + t519) * s12 - t668 * (-t952 - t121 + t917 + t433 - t87
     &6 - t954 + t953)
      t633 = intHLs16s25s26s34s56x1213D6eps1()
      t662 = s26 - s56 - s46
      t742 = -t662 + t75
      t800 = s26 + t199
      t801 = (s16 - s25 - s45 - t956 + t417 + t848 - t293) * s16 + s25 *
     & t742 - t379 * (-t1 * t800 + s34)
      t802 = t241 * s56
      t814 = t720 + t721
      t865 = t405 * t128
      t867 = -s45 + t405
      t904 = t351 * t405
      t921 = t444 * t1
      t928 = s34 - t417 - t177
      t971 = t134 + t287
      t565 = (-t16 * t715 + t26 * ((-s34 - t848) * s45 + t913) - t61 * t
     &814 + t722 + t723 + (t26 * (s34 - t848) - t61 * t417 - s26 + s45 -
     & t413 - s16) * s16 + s25 * (t932 - t956 + t568 - t200 - t142) + (-
     &t406 - t417 + t405) * s45 - t593 - t712 + t802 - t752) * s13 + ((t
     &175 - t925) * s16 - ((-ddreal(1) + t725) * s26 - t854 * s34 - t711
     & * (s45 - t200) + s25) * s25 - (t408 * s56 - t142 * t1 + t116 * t9
     &68 + t293 - t427 + t956) * s26 - t83 * (epinv * (s45 + t412 + t212
     &) - t212 + t293 - t412)) * s16 + (t134 - t287) * s26 + (epinv * t5
     &65 + t134 + t148 - t524) * s45 + s46 * t971 + s56 * t971 + t153 * 
     &t928 + t753 + t943 + t566 * t1
      t552 = ((-epinv * (s25 * (t199 + t522) - t121 + t148 - t366 + t770
     & + t554) - t121 + t148 + t770 - t870 - t442 + t751 * t552 - t177 *
     & s12) * s12 - t116 * (s26 * (t812 - t715) + t151 * (s26 + t848) + 
     &t208 * t595) + t15 * (t746 + t921) + t16 * ((t405 * t311 + t209 - 
     &t629) * s26 - t513 * t177) - t26 * ((s45 * t867 + t229 - t322) * s
     &25 + (-t720 + t133) * s46 + t66 * (t241 - t242 + t177) + t543 - t7
     &50 - t904 + t208 * (-t399 + t241)) + t61 * (s25 * ((s45 + s56 + t4
     &05) * s56 + t133) + (t716 + t865 + t316) * s46 + t148 * (s25 - t24
     &1) + t740 + t964) + t726 + t565 + ddreal(7) * t406 * t148) * s12 +
     & t668 * t801
      t554 = intHLs16s25s26s34s56x1231D6eps1()
      t565 = s34 - s45 + s46
      t593 = -t220 * t638 + t26 * (t391 + t325 + t879 + t87) - t877 + t8
     &81
      t595 = s34 * t800
      t638 = s34 - t417
      t723 = -s34 - t417
      t726 = t396 * s26
      t740 = t405 * s46
      t746 = t366 + t740
      t750 = epinv * t15
      t770 = t351 * s16
      t527 = -epinv * (s56 * (s46 * t698 - t524 * t16) + t116 * (s26 * (
     &t66 - t86) + (-t793 - t678) * s46) - t26 * (-t67 - t149 - t756 - t
     &692 + t795 - t444) + t61 * ((t121 + t527) * s56 + t687) - t609 + t
     &798 + s25 * (-t116 * t595 + t289 * t199 + t213 + t666 - t796) - (s
     &26 * t281 - t66) * s45 + t642 + t797) + (-t16 * t746 + t26 * (s45 
     &* t723 + t858) + t61 * (s46 * t638 - t397 * t127 + t228 * t84 + t2
     &30 + t736) - t897 - (-t406 - t724) * s25 - (t726 - t237) * s34 - t
     &134 - t712 + t220 * t788 - t513 * (t2 - t750) - t752) * s13 + s34 
     &* ((-s25 - s26 - s45 + s56 - s46) * s16 + t153 - t209 + t336 + t89
     &5)
      t724 = t147 * s16
      t353 = ((-epinv * (-t26 * (t86 - t570) + t630 + s25 * (t380 - t75)
     & + t121 - t655 + t85) + t565 * s56 + t85 - t870 - t90 + t751 * t62
     &6 - t896 - t268 * t48) * s12 + t116 * ((t892 + t316) * s46 + t770)
     & + t15 * ((t335 + t448) * s56 + t538) - t16 * ((t365 - t449) * s46
     & + s56 * (-s56 * t902 - t449) + t129 * t148 - t748 - t775) + t2 * 
     &((t366 + t85) * s56 + t724) + t26 * (s25 * (-s25 * t353 + t153) + 
     &(-s16 * t430 + t231) * s26 - t122 + t149 + t692 + t807 - t906) + t
     &61 * (s16 * t153 + (s34 * t563 + t449) * s45 + t444 + t392 * s16) 
     &+ t909 + t910 + t911 + t527) * s12 + t668 * t593
      t430 = intHLs16s25s26s34s56x1411D6eps1()
      t527 = t107 * t1
      t752 = t343 * epinv
      t800 = intHs16s25s26s34s56x1123D6eps1()
      t807 = intHL0s25s260s56x1022D4eps0()
      t881 = t379 * t320
      t892 = s25 + t84 + t43
      t897 = t320 * t375
      t371 = s12 * (t26 * s25 * t892 - s13 * t371 + t101 * s26 + t897) +
     & t12 * (-s25 * t742 + t323 + t881)
      t742 = intHL0s25s26s34s56x1120D2eps1()
      t902 = intHs16s25s26s34s56x1113D4eps1()
      t906 = t319 + t113
      t909 = s25 + s56 - s13 + s46
      t910 = s16 - s34 + s45
      t911 = t26 * t909
      t943 = (t910 + t911 + t289 + s12) * s12 + t668
      t964 = intHs16s25s26s34s56x1132D6eps1()
      t972 = intHL0s25s260s56x1021D2eps1()
      t973 = intHLs160000x0121D2eps0()
      t974 = intHLs160000x0122D4eps0()
      t975 = intHL0s25s26s34s56x1220D4eps0()
      t976 = s25 + s26 + s34 + s56
      t977 = intHL0s25s260s56x1012D2eps1()
      t978 = intHL0s25s260s56x1013D4eps1()
      t979 = intHL0s25s260s56x1010D0eps0()
      t980 = intHL0s25s260s56x1011D2eps0()
      t981 = intHL0s25s260s56x1020D2eps0()
      t982 = intHL0s25s26s34s56x1110D2eps0()
      t983 = intHs16s25s26s34s56x1113D4eps0()
      t984 = intHLs16s25s26s34s56x1311D4eps1()
      t985 = t223 + t225 + t75
      t986 = intHs160000x0111D0eps0()
      t987 = intHL0s25s26s34s56x1130D4eps0()
      t988 = intHs160s26s34s56x1022D4eps1()
      t989 = intHs16s25s26s34s56x1220D4eps1()
      t990 = intHs16s25s26s34s56x1321D6eps1()
      t991 = intHs16s25s26s34s56x1231D6eps1()
      t992 = intHs160000x0121D2eps1()
      t993 = intHs160000x0112D2eps1()
      t994 = intHLs160000x0112D2eps0()
      t995 = intHs16s25s26s34s56x1312D6eps1()
      t996 = intHs16s25s26s34s56x1222D6eps1()
      t997 = intHs16s25s26s34s56x1213D6eps1()
      t998 = intHLs160000x0122D4eps1()
      t999 = intHL0s25s260s56x1011D2eps1()
      t1000 = intHL0s25s26s34s56x1110D2eps1()
      t1001 = intHL0s25s260s56x1010D0eps1()
      t1002 = intHL0s25s260s56x1020D2eps1()
      t1003 = intHL0s25s260s56x1022D4eps1()
      t1004 = intHs16s25s26s34s56x1114D6eps0()
      t1005 = intHLs16s25s26s34s56x1411D6eps0()
      t1006 = intHLs16s25s26s34s56x1311D4eps0()
      t1007 = intHs160000x0111D0eps1()
      t1008 = intHL0s25s26s34s56x1310D4eps0()
      t1009 = intHLs160000x0112D2eps1()
      t1010 = intHs160s26s34s56x1031D4eps0()
      t1011 = ddreal(1) / s12
      t165 = ddreal(1) / t165
      t1012 = ddreal(1) / s25
      t1013 = ddreal(1) / t48
      t1014 = ddreal(1) / t25
      t1015 = ddreal(1) / t13
      t1016 = t1 * intHL0s25s26s34s56x1220D4eps1()
      t1017 = t1016 - t975
      t1018 = t1001 - t1002
      t1019 = intHL0s25s26s34s56x1120D2eps0() + intHL0s25s26s34s56x1210D
     &2eps0()
      t1020 = t999 + t1000
      t1021 = intHs160s26s34s56x1011D2eps0() + intHs16s25s26s34s56x1110D
     &2eps0()
      t1022 = ddreal(1) / (ecossin)**2
      t1023 = ddreal(1) / (gw)**2
      t1024 = t1 * t996
      t1025 = t327 * t8
      t1026 = t1015 * t10
      t1027 = t1018 * t1012
      t1028 = t1012 * t1020 * t27
      t1029 = t1012 * t1019
      t1030 = t1015 * t1011
      t1031 = t1023 * t1022 * propW34
      t1032 = intHs16s25s26s34s56x1120D2eps0() + intHs16s25s26s34s56x121
     &0D2eps0()
      t881 = t371 * t807 - t742 * (t405 * t189 + t46 * t48) + t972 * (-s
     &12 * (s13 * (s26 + s56 + t418) + s25 * (-t1 * t892 + t956) - t897)
     & + t12 * (s25 * (s26 - s45 + t932 - t241) + t323 + t881))
      t257 = -t370 * ((-t26 * (t208 + t376) + s25 * (s26 - t75) - t373 *
     & s26 + t374 + t383) * s12 - s13 * ((-t378 + s16) * s16 + t381 - t3
     &82) + t385) - t429 * t446 + (t443 * (-t26 * (t875 * t71 - t668 * (
     &s16 * (s45 + s56 + s46 - t241) + (s34 + t397 - t406) * s34 + t391 
     &+ t87 + t560 + t663)) + s12 * (s12 * (epinv * (t580 + t66 + t579 +
     & t133 + t587 + t584 + t583 - t586 + t585) - (epinv * t577 - s25 - 
     &s34 - s45 - s56 - t289 + t751) * s13 - s56 * t811 - t208 - t894 - 
     &t895 - t896) + t257)) + t561 * (((t116 * t177 * t531 - t26 * t812 
     &- (-t121 - t621 + t148 - t454 + t620) * epinv + s16 * (t1 * t375 +
     & s25 - ddreal(7) * t177) - s25 * t281 - t121 + t148 - t454 + t620 
     &- t621 - t870 - t617 + t310 * t177) * s12 + t644) * s12 + t668 * t
     &787 - t925 * t71) - t733 * t519) * t165 + t881 * t1013
      t374 = t53 * t55
      t375 = t327 * t464
      t376 = t18 * t78
      t381 = t8 * t1015
      t93 = t1015 * (t18 * (-t81 * (t26 * t92 - s13 * (t82 + t43) + s16 
     &* t83 + t85 + t86 + t87 - t88 + t93) + t1011 * t257) + t21 * (-t37
     &5 + t374) + (t342 * (s12 * (s46 * t202 + s56 * t281 + t345 * s56 +
     & t321 * t58 + t148 + t344 + t66) + t346) + t347 * (-t361 * s34 + s
     &12 * (-(t360 + t359) * s34 - t16 * t356 + t26 * t357 - s13 * t355 
     &+ s16 * (t175 - t293) + t351 + t358) + t12 * t354) + t395 * ((t26 
     &* ((s26 + t397) * s46 + t325 + t391) + t93 - s13 * (s16 + s26 + t3
     &97 + t398 - t241 + t43) + s16 * t339 + (s26 - t241) * s25 + (t228 
     &+ t400 + t399 - t241) * s46 + (t228 + t399 - t241) * s56 + t230 + 
     &t401 + t402 - t403 - t341) * s12 + t80 * t25) + t404 * t421) * t10
     &11 * t8 + (-t468 * t474 - t475 * t482) * t1014 * t10) * t1012
      t130 = -t508 * (((t26 * (t514 + t513 + t210) + (-t511 - t512 + s13
     &) * s13 + t121 + t148 - t204 + t66 + t516 + t517 + t518) * s12 + t
     &116 * s46 * t541 + t15 * t366 * t127 - t16 * t544 + t26 * (t134 * 
     &t537 - t153 * t373 + t456 * t66 + t538 - t539 + t540) + t61 * ((t8
     &6 + t534) * s25 + (-t66 + t133 - t436) * s26 + s46 * t532 - t535 -
     & t536) - (-t26 * (t391 + t134) + t61 * (-t123 + t523 - t338) + s16
     & * (-t520 + t139 + t521 - t289) - s25 * (t199 - t522) - s45 * t343
     & - t148 + t524 - t525 - t549) * s13 + s16 * t530 - s34 * t546 - t1
     &22 - t444 - t550 - t335 * t548 - t316 * t128) * s12 + t196 * t509 
     &* t510) + t558 * t562 + t660 * t130
      t139 = t611 * t669 + t845 * ((t699 * s12 + (s46 * (t874 + t803) + 
     &ddreal(11) * t366) * s56 + t116 * ((-t241 * s26 + t716) * s46 + t3
     &16 * t531 - t86 * t194) + t15 * s26 * ((s25 + s26) * s46 + ((s34 -
     & s45) * epinv + s25) * s56 + t335) - t16 * (epinv * t544 + t207 * 
     &t456 - t67 + t776 + t872) + t26 * (((-s34 + t417) * s46 - t415) * 
     &s45 + (epinv * (t335 - t134) + t664) * s46 + t153 * (-epinv * t373
     & + s25) + t66 * (s25 + t708) - t873 - t415 * t246) + t61 * (((s34 
     &+ t860) * s34 - t861) * s25 + (s56 * (-t228 + t242) - t133 - t852)
     & * s56 - t569 + t793 * t228 + t847 * t532) - t528 + ddreal(7) * s5
     &6 * t871 + ddreal(7) * t208 * t143 + ddreal(7) * t693 + ddreal(10)
     & * s56 * (t208 + t513) - t645 * s26) * s12 + t72 * t11 * t509 * t5
     &10)
      t175 = t186 * t215 - t222 * (t26 * ((s25 + s34 + s45 + t228) * s45
     & + (-s26 + s34 + s45 + t177) * s46 + (s12 - s26 + s34 + s45 + t177
     &) * s56 + t229 + t230 + t231 + t232) + t61 * epinv * (s45 * t233 +
     & t234) + epinv * (t205 + t151 + t121 + t206 + t153 + t133 + t86) -
     & s13 * t226 + t235 * t175 + t213 * epinv)
      t186 = -t186 * t214 + t216 * (-s12 * (t26 * (s25 + s46 - t217) + t
     &142 + epinv * t218 + s26 + s34 + t219 - t220 + t221) + t1 * (t70 +
     & t196))
      t229 = t244 * t21 * (t264 + t265 + t266 - t268 * (s26 * t98 - (t13
     &2 + t60) * s45 + (-s25 - s16 - s45 + s13 - t60 - s46) * s46 + (-s2
     &5 + s12 - s16 - t63 - t62) * s56 + t258 - t259 + t260 - t68) + t26
     &7)
      t17 = epinv2 * (t3 * (t5 * (-t16 * ((-s25 * t836 - s26 * (-s13 + t
     &113 + t837) + t258) * s12 + (-s12 * (s13 - t256) + t691 - t276) * 
     &s34 + (-s12 * t839 + t821) * s45 + ((t279 + t296 - t840 + s56) * s
     &56 - (s13 + t840) * s16 - (-s16 + t106) * s25 + (t279 - t841 + s26
     &) * s26 + (t842 - s34) * s34 + (s25 + s12 - s34 - s13 + t102 + t28
     &9 + s45) * s45 + t162 + t830 - t361) * s56 + t831 + t832 + t833) -
     & ddreal(12) * t835 + t17 * t25 * (s12 * t245 + (-t822 - s45) * s45
     & + (-s46 - t823) * s46 + (t26 * t843 - s56 - t820) * s56 - t276 - 
     &t834 * t187)) + t19 * ((-s13 * t25 + s56 * t676 + t326 + t453 + t6
     &6 - t661 - t678 - t88) * s12 - t346)) + t844) * t1013 * t1014
      t257 = (intHs160s26s34s56x1011D2eps1() + intHs16s25s26s34s56x1110D
     &2eps1()) * t27
      t276 = intHs160s26s34s56x1020D2eps1() * t1
      t279 = t1031 * t1015
      t321 = t279 * t1012
      t17 = t321 * (t1011 * (t12 * (t11 * (t10 * t9 + t14 * t21 + t30 * 
     &t8) + t22 * t10 * t23) + t165 * (t10 * t139 + t18 * (-t575 * (t26 
     &* (-t328 * t510 + t48 * t71) - s12 * (s12 * ((-t577 + s13) * s13 +
     & t133 + t66 + t579 + t580 + t587 + t584 + t583 - t586 + t585) - t1
     &16 * s45 * (t231 - t121 + t602) + t16 * (t614 - t612 + t613 + t444
     &) + t2 * t600 - t26 * (s45 * (t66 + t86) + (s45 * t603 - t66) * s4
     &6 - t148 * t300 + t604 - t605) + t61 * (-t134 * t246 + t569 + t610
     &) + (t599 + t549) * s13 + (-t211 + t210) * s25 + (t606 - t153) * s
     &34 + t550 - t607 + t608 + t615 + t428 * s26)) - t616 * t542 + t799
     & * t91) + t8 * t130) + t18 * t175 + t8 * t186 - t229 - t17) + t101
     &4 * (t10 * ((-t28 - t29) * s34 + t257) - t276 * t8) * t24)
      t91 = t188 + t191
      t130 = t1 * t10
      t139 = t35 * (s16 + s26 - s34 - s45 + t36)
      t175 = t10 * t38
      t186 = t241 * intHL0s25s26s34s56x1210D2eps1()
      t229 = s34 * t24
      t328 = t10 * (t1 * (t229 * t37 + t482 * t989) - t166 * ((epinv * t
     &24 + s34 - t168) * s34 - t170 + t171 + t172 + t173)) + t8 * (t174 
     &* (t177 * t24 + t176) + t24 * t96)
      t339 = intHs16s25s26s34s56x1210D2eps1() * t1
      t346 = t1011 * t18
      t205 = t279 * (t1012 * (t1011 * (t12 * ((t139 * t21 + t175 * t40) 
     &* s34 + t11 * (epinv * t8 * t31 + t21 * t33 * t34 - t130 * t32)) +
     & t18 * (-t197 * (-s13 * t201 + t26 * (t208 + t209 + t66 + t148 + t
     &207) + t61 * (t211 + t210) + t121 + t133 + t151 + t153 + t204 + t2
     &05 + t206 + t86 + t213) + t236 * ((t26 * (s25 - s34 + s56 + s46 - 
     &t237) - (s25 + s34 + s56 + s46) * epinv + s16 * t34 + t219 + t238 
     &+ t239) * s12 - t12 * t226) - t186 * t47 * t1013) + t8 * (t178 * t
     &179 + t180 * t182 - t192 * ((s12 - t193 - t195) * s12 + t196) + t2
     &40 * ((t36 - t26 * (s25 + s56 + s46 - t237) - t142 + (s25 + s56 + 
     &s46) * epinv - s26 - t219 + t220 - t241) * s12 + t12 * (t223 * t26
     & + t225 - t241 - t242 - t243)))) + t1014 * t328) - t339 * t8 * t24
     & * t1014 - t346 * t1013 * t91 * t189)
      t206 = epinv * intHLs160000x0121D2eps1()
      t213 = t27 * intHLs160000x0111D2eps1()
      t219 = s16 * (t1 * (-t998 + t1009) - t973 + t974 - t994 - t206) + 
     &t183 + t213
      t226 = t1 * intHs16s25s26s34s56x1114D6eps1()
      t328 = t226 - t1004
      t358 = t1 * intHs16s25s26s34s56x1141D6eps1()
      t383 = t430 * t18
      t385 = t1012 * (t358 - t56)
      t401 = t1012 * t1011
      t403 = t1015 * t1022 * propW34
      t23 = t403 * (t1023 * (t1011 * (t21 * (t1012 * (-t249 * t483 + t26
     &9 * (-epinv * t249 + t264 + t265 + t266 + t267) + t327 * t23 * t10
     &05) - t196 * t466) - t383 * s34 * ((s12 + s16 - s34 - t752 + t418 
     &- t751 + t417 + t848 + t60) * s12 - t12 * (s34 + s45 + t223 - t527
     &)) * t1012 + t1025 * t23 * t1012 * t328) + t385 * (s26 * t57 + (t5
     &9 + t60) * s45 - (-s25 - s16 + s34 - s45 + s13 - t60 - s46) * s46 
     &+ (s25 - s12 + s16 - s34 + t63 + t62) * s56 + t64 - t65 + t68) * t
     &10) + t401 * t219 * t187 * t185 * propW16)
      t63 = t1 * intHL0s25s26s34s56x1130D4eps1()
      t64 = -t63 + t987
      t68 = t1 * t990 - t486
      t92 = t8 * t422 * (((-t405 + t423) * s12 - t16 * t92 + t26 * (s26 
     &* t426 + (s45 + s56 + t418 + t417) * s34 - t207 - t390) - t61 * t4
     &05 * t320 + epinv * (t357 + t121 + t367 + t66 + t338) - s13 * (t26
     & * (-s46 + t241) - t221 + t424) + s16 * (epinv * (s26 - s34 + s45 
     &+ s46 + s56) + t413 - t425) - t428) * s12 + t369 * t72)
      t249 = t1 * t978
      t265 = t1 * t461
      t266 = t1 * t165
      t267 = t983 * s34
      t196 = t1026 * (t1012 * (-t267 * t943 + t266 * (t564 * t997 + t573
     & * t995)) - t196 * t465)
      t45 = t1012 * (t1015 * (t12 * (-s34 * t906 * t1006 + t45 * t64) + 
     &t272 * t485 - t274 * t484) + t265 * t12) * t21
      t272 = t1026 * t1012
      t327 = t272 * t165
      t108 = t327 * t68 * (t26 * t67 - (-t99 * s25 - (-t283 + t111) * s2
     &6 - t113 * t488 + t489) * s26 - (-s13 * t270 + (t473 + t113) * s25
     & - t231) * s34 - ((s13 + t490) * s25 + (s25 - t491) * s26 + (t108 
     &- t111) * s34 - t119 - t492 - t507) * s45 - ((-t26 * t290 - s46 - 
     &t289 - t498) * s46 - (-t493 + t494) * s26 - (t495 + t118 + t102 - 
     &t104 + s45) * s45 - t119 + t259 - t496 - t497) * s46 - (-(t290 * t
     &61 + t102 + t495 + t505) * s56 + t503 * s13 + (s12 - t499 + t104) 
     &* s26 + (-t495 - t360 - t500 - s45) * s45 - (t116 * t290 + t109 + 
     &t501 + t502) * s46 + t259 - t496 - t504 - t492) * s56 - t100 * t50
     &6 + t299 * t487)
      t270 = t1026 * (t1012 * (t163 * t991 * t165 + t12 * t463) + t24 * 
     &t460 * t1014) * t1
      t421 = t1 * t800
      t428 = -t421 + t330
      t126 = t165 * (-t564 * t571 - t622 * (((t88 - t26 * (t259 - t570 +
     & t316) + (t626 - s13) * s13 + t148 + t334 + t581 + t66 - t629 + t6
     &30) * s12 - t116 * ((t793 + t231) * s46 + t600 - t613 + t794 + t23
     &1 * s56) + t16 * s56 * (t209 - t524) + t26 * (t128 * t148 + t149 +
     & t444 + t67 + t695 + t739 + t756 - t795) + t61 * (t335 * t199 + t6
     &12 + t748 + t778) - (-s16 * t788 - (t161 - t526) * s46 + t26 * (t1
     &47 + t366) - t61 * (t790 + t121 - t66 - t148 - t90) - t555 + t789 
     &- t86 + t549) * s13 - s25 * t796 + t536 - t775 + t797 + t86 * t127
     & + t798 - t609 + t513 * t698) * s12 - t12 * ((t785 * s16 - t16 * (
     &s26 * t83 + t513) - t26 * (-t153 + t66 + t148) - t357 + ddreal(9) 
     &* t650 + t656) * s16 + t767)) - t700 * t126 + t926 * (s12 * t707 +
     & t668 * t701) + t633 * t552 + t554 * t353) + t318 * t329
      t353 = s34 * t902
      t337 = t381 * (t428 * ((-t26 * (s26 * t337 - t338) - t61 * t340 - 
     &s13 * (t319 - t144) + s16 * t333 + t148 + t66 + t334 + t335 + t336
     & - t341) * s12 + t332 * t319 * t12) - t362 * ((-t368 - t26 * (-t36
     &6 + t356) - t365 * t61 - s13 * t364 + s16 * (s25 + s26 - s34 + s45
     & + s46 + s56) + t121 + t338 + t357 + t66 + t367) * s12 + t369 * t1
     &2) - t353 * (-epinv * t943 + t12 * t906))
      t356 = t1 * t964
      t251 = t401 * (t21 * (-t1015 * (t278 * (epinv * t274 + s25 * t305 
     &+ (-t61 * s25 * (s25 + t281) + s13 * t280 - t282 - t277) * s45 + (
     &-t26 * (-t253 + t121 + t287) - (s16 + t284) * s25 - (t251 * t26 + 
     &t101) * s26 - (s25 - t285) * s34 + (-t26 * t286 - s13 + s45 - t105
     & - t109) * s45 + t306) * s46 + ((-t251 - t296 - t101 - s56) * s56 
     &- t26 * s25 * (s25 + t283) + (s13 - t110 - t102 - t289 - s45) * s4
     &5 + (-t290 - t291 + t293 - t238 - t292 - s46) * s46 + t253 - t254 
     &- t294) * s56 - t264 * t61) + t307 * t317) - t12 * t51) + t337)
      t126 = t1031 * (t21 * t95 + t1026 * (t1012 * (-t12 * t52 + (-s12 *
     & (t16 * (t325 + t392) + t26 * (t389 + t390) - t61 * (-t391 + t340)
     & - s13 * (-t387 * t61 + t388) + s16 * (t26 * t333 + s56) - t86 + t
     &393 * s56) - t394 + t361 * t320) * t1011 * (-t356 + t386)) - t69 *
     & (t11 * (s12 + t72) + t70) * t1011 - t50 * t24 * t1014) + t251 + t
     &1030 * (t73 * (t12 * t76 + t77) + t1012 * t126) * t18)
      t251 = t1 * intHs160s26s34s56x1013D4eps1() - intHs160s26s34s56x101
     &3D4eps0()
      t254 = t1030 * t18
      t264 = intHs160s26s34s56x1031D4eps1() * t1
      t274 = t1015 * t1012
      t88 = t274 * (t10 * (-t1011 * (t551 * t573 + t564 * t768) - t163 *
     & t97) + t346 * (t559 * (((-t88 + t26 * (t524 + t316) + t61 * (t325
     & + t259 - t323) + (-t601 + s13) * s13 - t121 + t590 + t625) * s12 
     &+ t116 * s34 * ((s25 + s46) * s46 + t207) + t26 * (s25 * t755 - (t
     &134 + t66 + t148) * s26 + (t655 - t148 + t86) * s45 - t756) - t61 
     &* (t686 + t748 + t739 + t754) + ddreal(7) * t600 - (t16 * (t151 - 
     &t338) + t61 * (t365 + t121 - t596) + s16 * t653 - t134 - t153 - t3
     &92 + t555 - t666 - t549) * s13 - (t457 + t757 + t153 + t66 + t148)
     & * s25 - t149 - t607 - t67 - t689 - t690 - t758 - t759 - t760 + t7
     &61 - t762 + t636 * t524) * s12 + t12 * (s16 * (-s16 * t26 * t763 -
     & t765) + t767)) - t573 * t574)) * t165
      t88 = t1031 * (t88 - t254 * (t241 * t47 * t49 * t1012 + t189 * t19
     &0) * t1013 + t264 * t8 * t1012 + t272 * t229 * t251 * t1014)
      t45 = t88 + t126 + t1031 * (t1011 * (t1015 * (t1012 * (t18 * ((t98
     &4 * ((s12 + t26 * (-epinv * t440 + s25 + t417 - t751 + t848) + t28
     &9 + t910) * s12 - t12 * t985) + t47 * t1008 * t1013) * s34 - t452 
     &* (-epinv * t329 + t15 * t366 * t84 + t16 * ((t147 + t338) * s16 +
     & s26 * t455 + s56 * (s46 * t434 + t351)) - t26 * (s13 * ((s12 + s2
     &6 - s34 + s56) * s46 + (s12 + s26 + s56) * s56 + t325 + t87 + t453
     &) + (-t26 * t455 - t61 * t338 - s12 * t199 - s26 * (s34 + t425) - 
     &t148 - t153 - t453 - t454 - t87) * s25 - t121 * t194 + t148 * t459
     & + (-t366 - t457 - t66) * s12 + (-t366 - t208 - t392 - t153 - t66)
     & * s16 + (t151 - t153) * s26 + (-t153 + t458) * s56 - t444 - t67 +
     & t85 * s34))) - t92) + t249 * t18 * t189 * t1013) + t196 + t45) - 
     &t1010 * t8 * t1012 + t108 + t270)
      t88 = t1 * intHLs160000x0113D4eps1() - intHLs160000x0113D4eps0()
      t92 = t1030 * propW16
      t108 = t92 * propW34 * s16 * t185
      t126 = -ddreal(2) / ddreal(3)
      t163 = ddreal(2) / ddreal(9)
      t196 = ddreal(1) / ddreal(9)
      t229 = s25 - s26 + s56 - s13 + s46 - t275 + t142
      t270 = intHs16s25s26s34s56x1411D6eps0()
      t280 = t101 + t60
      t286 = s25 - s12 - s16 + s26 + s34 - s13 + s46 + t54
      t290 = (t102 + t106) * s25
      t296 = (s12 + s13 + s16 + s25) * s26
      t305 = (s26 + t101) * s34
      t317 = t280 * s45
      t323 = (-s12 - s16 + s25 + s34 + s45) * s46
      t329 = (t286 + s56) * s56
      t337 = (s16 + t104) * s25
      t340 = s25 + s16 - s45
      t346 = s12 + s26
      t365 = t26 * t346
      t367 = s25 - s12 - s16 - s13 + s46 + t54 + t144
      t390 = (t58 + t490 + t289) * s34
      t394 = t16 * t134
      t425 = t26 * (t170 + t390)
      t434 = s26 * t132
      t440 = t132 * s16
      t446 = t113 * s16
      t455 = s26 + s45 - s13 + s46
      t457 = s12 + s13 + s16
      t482 = t128 - t275
      t511 = s25 + s26 - s13 + s46
      t514 = t511 - t275 + t142
      t518 = t280 * s34
      t519 = t60 * t457
      t523 = t16 * s13
      t525 = s12 + s16 + s26 - s34
      t528 = t26 * t525
      t530 = s16 + t113
      t532 = -s16 + s26 + s34 + s45 + t111
      t511 = t26 * t511 + t142 - t273
      t535 = (s25 + s12 + s13 + t102) * s26
      t540 = (s16 + t105 + t60) * s34
      t542 = (s16 + t105 + t289) * s45
      t544 = t530 * s25
      t546 = t26 * (t66 - t544)
      t549 = (s16 + t523) * s25
      t550 = (t107 + t105) * s45
      t552 = (-s45 + t528 + t101) * s46
      t555 = t26 * t673
      t562 = s45 + s13
      t564 = s34 * t387
      t573 = s45 + s56 - s13 + s46
      t577 = s25 + s26 + s46 + s56
      t579 = s12 * t577
      t580 = s13 * s26
      t583 = (s56 + t455 + t111) * s56
      t584 = s25 - s16 + s56 - s13
      t585 = s12 - s34 - s46
      t586 = -t26 * t585 + t142 + t584
      t596 = t61 * t273
      t585 = t26 * t584 + t142 - t585
      t599 = -s12 - s16 + s26 + s34 + t911 + t142
      t600 = t61 * t58
      t601 = s25 + s12
      t602 = s12 + s16 + s25
      t604 = t26 * t134
      t607 = s25 - s16 + s26
      t608 = t457 * s26
      t273 = t273 * s46
      t609 = s25 + s26 + s34
      t610 = (s16 + t478) * s34
      t615 = t609 * s46
      t286 = (t286 + s56) * s56
      t617 = t478 * s45
      t621 = s12 + s25 + s26
      t626 = t74 * t246
      t629 = (s56 + t42) * s56
      t630 = -s12 * t115 + s46 * t621 + t626 + t629
      t644 = s34 * (t459 + s34)
      t525 = -s46 * t525 - t170 - t629 + t644
      t629 = (-s26 - t600) * s34 + (-t297 - t144) * s45 + (t340 - t144 +
     & t365) * s46 + (s56 + t42 - t144) * s56 + t434 + t440 + t446
      t645 = s12 + s25 + s26 - s34
      t653 = s12 * s16
      t656 = t673 * s45
      t666 = t645 * s46
      t668 = s25 * t107
      t691 = (-s13 + s25 + s26 - s34 + s45 + s46 + s56) * s56
      t695 = t26 * t12
      t698 = t26 * t11
      t699 = s12 - s34 - s45 - s46
      t707 = s25 - s45 + s56 - s13 + s46 + t275 + t289
      t755 = s25 + s26 + s45 + s56 - s13 + s46 - t39
      t756 = epinv * t699
      t758 = (t602 + s26) * s26
      t372 = t74 * t372
      t530 = (-s46 + t530) * s25
      t760 = t530 + t758 + t171 - t372
      t470 = t26 * t644 + t151 + t171 + t172 + t173 - t470
      t644 = t405 * t699
      t761 = t58 * s25
      t762 = -s13 * t48 + (t84 + t423 + s25) * s25 + s26 * t295 - s34 * 
     &t843 + t134 + t338 + t66 + t366 * t26
      t767 = t117 - t109
      t500 = t500 * s16
      t778 = s13 * t132
      t788 = (-t105 - t114) * s25 - t61 * (t123 - t778) + s12 * (s13 - t
     &837) + t506
      t790 = s25 * t246
      t796 = s13 * (s16 + t111)
      t797 = t16 * t123
      t798 = t113 + t99 + t289
      t820 = s25 - s12 - s16
      t821 = s25 * (s25 - t457)
      t822 = (t840 + t104) * s16
      t823 = (s13 + t111 + t159 + t106) * s26
      t830 = t116 * t123
      t831 = ddreal(7) * t107
      t832 = -t116 * t704
      t827 = t26 * (s25 * t827 + t70)
      t833 = t116 * t311
      t835 = t61 * (s16 + s45 - s13)
      t843 = s13 * (s25 - s34 - s45 + s46 + t164)
      t844 = t107 * s45
      t860 = t15 * (t209 + t121)
      t870 = ddreal(9) * s26 * t146 + ddreal(9) * t449
      t872 = t116 * s16
      t881 = s13 + t159
      t892 = ddreal(9) * s12
      t895 = t61 * (t134 + t299)
      t896 = (t872 + t841) * s16
      t897 = t2 * s12
      t906 = t15 * t107
      t910 = t116 * t100
      t932 = (s26 - t54 + t596) * s56 + (t892 + t159) * s16 + t26 * (t12
     &1 + t153 + t70) + t895 + (s25 + t117 + t897 + t156) * s26 + (s13 -
     & t101 - ddreal(10) * t107 - t826) * s34 - (t906 + t111 - t161 + t1
     &06) * s45 + (s26 - t54 + t910 - t670) * s46 - t162
      t943 = t258 * t824
      t1033 = ((t881 + t892) * s16 + t361) * s25
      t1034 = t104 * t123
      t895 = ((s25 + t523 + t826 + t114 + s26) * s26 + s12 * (s13 + t291
     & + t908) + s16 * (t159 + t103) - t489 + t361) * s26 + ((t476 + t25
     &6) * s34 + s12 * (-t476 - t576) + (-t102 - t672 - t631 + t104) * s
     &16 - s26 * (t110 + t138)) * s34 + ((s34 + t478 - t109) * s45 - t16
     & * (t121 + t479) + t830 - (s13 + t105 + t292) * s26 + (-t297 * t61
     & + s34 + t291 - t98) * s34 + t822) * s45 + ((t490 - t144) * s46 - 
     &(-s25 - t117 - t841 - t114 - s26) * s26 + (s13 - t826 - t697 - t10
     &1 - t908) * s34 - (-t26 * t227 + t100) * s45 - t162 + t896 + t361 
     &+ t895) * s46 + t932 * s56 - t943 + t1033 - t1034
      t932 = (-s12 - t291 - t99 - s26) * s26
      t1035 = (-t798 + s34) * s34
      t1036 = (-t137 * t26 + t16 * t311 + t106 + t138 + t872) * s56
      t255 = t255 * t148
      t1037 = t167 * t287
      t1038 = t506 - t70
      t1039 = t113 * s13
      t1040 = t1038 + t1039
      t1041 = s13 - t106
      t1042 = ddreal(7) * t123
      t1043 = t61 * t601
      t1044 = ddreal(13) * s16
      t1045 = t26 * t248
      t1046 = ddreal(7) * t100
      t1047 = t16 * t673
      t1048 = t287 + t70
      t1049 = (s13 + t841) * s16
      t1050 = (-t117 + t99 + t106) * s25
      t1051 = (t113 + t1047 + t872 - s45) * s45
      t1052 = t121 + t134
      t1053 = t61 * t1052
      t1054 = t26 * t601
      t1055 = t11 + t1054
      t1056 = t15 * s12
      t1057 = s13 + t1056
      t1058 = t1057 * s16
      t1059 = (s12 - t161 + t872 + t289 - s45) * s45
      t1060 = t61 * t134
      t1061 = s13 + t872
      t1062 = (s13 + t897) * s16
      t1063 = (s13 + t1056 + t114) * s25
      t1064 = ddreal(10) * s25
      t1065 = t61 * (t123 - t580)
      t1066 = s16 + t285 + t105
      t1067 = t476 * s16
      t1068 = (t840 + t159) * s25
      t457 = t26 * (-t121 + t134 - t1067) - (s12 + t145 + t291) * s26 - 
     &(t26 * t457 + s26 - t110) * s34 - t1068
      t1069 = s25 + s12 - s34
      t1070 = t331 + t113 - t144
      t1071 = t16 * t601
      t1072 = (t145 + t1043 + t117 + s26) * s26
      t1073 = (s26 + t1045) * s45
      t1074 = t1070 * t66
      t1075 = t1070 * t148
      t1070 = (t26 * (s46 * t1070 + t121) + t1060 + t797 + (s13 + t145 +
     & t111 + t826) * s26 - (-t117 + t1071 + t805 + t114) * s34 - (s12 +
     & t502 + t872 + t289 - s45) * s45 + t1058 + t290 + t828) * s56
      t1076 = ((t1061 + t897) * s16 + t361) * s25
      t1077 = ((s13 + t111 + t1046) * s26 - (-t117 + t105 + t826 + t697 
     &+ t114) * s34 - (t1069 * t26 - s45 + t238 + t872) * s45 + t1049 - 
     &t162 + t290 + t797 + t361 + t1053) * s46
      t1078 = (-t506 + t70 + t1063 + t1062 + t1042 + t1072) * s26
      t1079 = (-t1041 * s16 - t1040) * s16 + (t1066 * s34 + s12 * (-t110
     & - t491 - t289) + (t523 - t1064) * s16 - t684 - ddreal(9) * s26 * 
     &t132 - t1065) * s34 + (t457 + t1073) * s45 + t1074 + t1075 + t1070
     & + t1076 + t1077 + t1078 + t783
      t1080 = (t360 + t840 + t159) * s34
      t1081 = (-t1041 * s16 - t1040) * s16
      t1082 = t16 * t287
      t1083 = (s13 + t826) * s16
      t1084 = t26 * t287
      t1085 = s25 - s12
      t1086 = t283 * s16
      t1087 = t16 * t100
      t824 = (t824 + s16) * s16
      t1088 = (-t117 + t292) * s25
      t1089 = t379 + t1054
      t1090 = t117 + t826
      t1091 = t16 * t107
      t1092 = s13 + t113
      t1093 = t61 * t123
      t1094 = -s13 + t291 + t60
      t1095 = s13 * (s26 + t102)
      t1096 = t114 * s25
      t1097 = (t99 + t618 + t289) * s26
      t252 = -t252 - t101
      t331 = t331 + t298
      t1098 = (t1043 + t117 + t99 + s26) * s26
      t1099 = t1092 * s16
      t1100 = t331 * t66
      t1101 = t331 * t148
      t1069 = (t26 * (-(t1069 + t164) * s45 + t134 + t70) - t61 * t704 +
     & (s13 + t145 + t111 + t841) * s26 + (s13 - t840 - t138 - t872 - t1
     &01) * s34 + t1058 + t153 - t162 + t290) * s46
      t331 = (t26 * (s46 * t331 + t121 + t134) + t1093 + (s13 + t910 + t
     &111) * s26 + (s13 - t360 - t1043 - t872) * s34 - (-s45 + s12 + t99
     & + t502 + t289) * s45 + t1083 + t290 + t828) * s56
      t1102 = (-t123 - t1099) * s13
      t1063 = (-t506 + t70 + t1063 + t1083 + t1093 + t1098) * s26
      t1083 = -(s16 + t117) * s25 + (t488 + t840 + t114) * s26 + t70 + t
     &1083 + t830 + t1084 + t492 + t1035
      t896 = t61 * s13 * t115 - t26 * s25 * t100 + t1038 - t1082 - t1098
     & - t896
      t169 = (t169 - t111) * s45
      t487 = -t169 + t26 * ((-t487 + s25) * s25 + (t1085 - t109 - t60 + 
     &s34) * s34 + t121) + t797 + (s12 + t159 - t101) * s26 + t1086
      t1098 = t124 + t1037
      t119 = t26 * t1098 - epinv * ((-s12 * t1094 - s34 * t252 + t1095 -
     & t1096 - t1097) * s34 + (-s45 * t839 - t26 * (-t120 + t121) - (s12
     & + t291 + t109) * s26 + (s34 - t98 + t291 - t102) * s34 - t1068 + 
     &t123) * s45 + t1076 + t1100 + t1101 + t1069 + t331 + t1102 + t1063
     &) - (-t70 + t119 - t829) * s16 - (s13 * t836 - s16 * t493) * s25 -
     & t896 * s26 - t1083 * s34 - t487 * s45 + (t26 * (t824 + t134 + t28
     &7 + t70) + t492 + (s13 + t110 + t99 + t841) * s26 - (t118 + t111 +
     & t1087) * s34 - (t113 + t238 + t872 - t144 - s45) * s45 - t162 + t
     &1088) * s46 + (t26 * (s46 * t1089 + t121 + t123 + t134 + t287) + s
     &16 * t1090 + (s13 + t110 + t826 + t109) * s26 - (t111 + t1091 + t1
     &06) * s34 - (t61 * t673 + s12 - s45 + t872) * s45 + t1088 + t828) 
     &* s56 + t148 * t1089 + t66 * t1089
      t487 = s12 - t117
      t836 = t487 * s16
      t896 = t104 * s12
      t1083 = -t1038 - t896
      t1088 = ddreal(9) * t123
      t1089 = s12 * (t107 * t2 + t58)
      t1103 = ddreal(16) * s26
      t1104 = -s25 + s26
      t1105 = ddreal(10) * s16
      t1106 = ddreal(9) * t107
      t1107 = (s25 - s13 + t113 + t109) * s25
      t1108 = (t145 + t897) * s16
      t275 = s25 + s26 - s45 + t275
      t1109 = (t145 + t841) * s16
      t1110 = (t58 + t840 + t831) * s34
      t1111 = (s25 + t117 + t114 + t106 + s26) * s26
      t1112 = t1083 * s16
      t1113 = (-t100 * (s13 - t109) + t479) * s25
      t1114 = t283 * t61
      t275 = -t26 * (-t124 + t135) + ((t111 + t99) * s25 - (-s12 - s25 -
     & t908) * s12 - (s13 + t132) * s13 + t1088 + t1111) * s26 + (-(t110
     &3 + t101) * s16 - t26 * (t335 + t70) + ddreal(7) * t704 + s13 * (t
     &246 + t109) - t287 - t1089 + t1110) * s34 + (-(s16 + s25 - s26 + s
     &34) * s45 + (t360 - t144) * s34 + t116 * s16 * t531 - t26 * s26 * 
     &t1104 + s12 * (s25 + s16 - s26 + t144) + s13 * (t495 + t109) - t28
     &7) * s45 + ((s13 + t841 + t1105 + t101) * s26 + (s13 - t1056 - t11
     &1 - t1106) * s34 - t1051 - t162 + t1107 + t1108 + t361 + t394 + t4
     &92) * s46 + (t26 * (t275 * s46 + t121) + t394 + (s13 + t826 + t101
     & + t908) * s26 + (s13 - t111 - t826 - t501 - t908) * s34 - t1059 +
     & t1107 + t828 + t1109) * s56 + t148 * t275 + t66 * t275 + t1112 + 
     &t1113 + t1114 * t123
      t1115 = s13 - t113
      t1116 = ddreal(10) * s12
      t1117 = (-t115 + s13) * s13
      t1118 = ddreal(12) * s25
      t1119 = s12 + t117
      t1120 = (s16 + t1119) * s16
      t1121 = s25 - s26 + s34
      t1122 = ddreal(7) * t86
      t1123 = t26 * (t134 + t70)
      t1124 = t16 * t283
      t1125 = t16 * t476
      t1126 = s34 - s13
      t1127 = t26 * t1126
      t1128 = (t826 + t114) * s25
      t1129 = (t132 * t16 - t117 - t670 + t805 + t826) * s46
      t1130 = t61 * (-t1086 + t162)
      t1131 = t26 * (t134 - t153 + t148)
      t1132 = s25 + t109
      t1133 = (t1132 + s13) * s13
      t1134 = ddreal(16) * s16
      t1135 = -s12 * (s13 - t99) - (s13 - t872) * s16 + t479
      t528 = t1135 * s25 + ((t291 + t840 + t1105) * s26 + (t111 + t1105)
     & * s25 + t1088 - (-s12 + s13 - t101 - t908) * s12 - t1133) * s26 +
     & (-(t501 + t1134) * s26 - t16 * t668 - t1042 - t361 + s13 * (s25 +
     & t109 + t60) - t1089 - t287 + t1110) * s34 + (-s45 * t1121 + t61 *
     & (-t121 - t134 + t258) - t125 - (s25 + t503) * s25 + (-t283 - t159
     &) * s26 + (-s13 + t906 + t106) * s34) * s45 + (t26 * (t172 + t121 
     &+ t134) + (s25 - s13 + t145 + t106) * s25 + (t113 + t99 + t618) * 
     &s26 + (s13 - t113 - t238 - t837) * s34 + (-s25 + s12 - s26 + s34 -
     & t102) * s45 - t162) * s46 + ((t592 + t528) * s56 + (s25 - s13 + t
     &113 + t99) * s25 + (t1056 + t291 + t156) * s26 + (s13 - t111 - t82
     &6 - t1106) * s34 + (-s12 - t145 - t1047 + s45) * s45 + t552 + t828
     & + t1109 + t394 + t684) * s56 + t1112
      t447 = -t26 * (-t124 + t135 - t447) + t61 * (t123 * t283 + t122) +
     & t528
      t528 = (t16 * (s16 - s34 + s46) + t106 + t138 + t494) * s56
      t645 = t645 * t148
      t1047 = -t26 * (-t112 * t134 + t645 + t67) + epinv * t447 - (-s16 
     &* t1115 - t1040) * s16 - ((s13 + t908 + t1116) * s16 + t361) * s25
     & - ((s12 + t102 + t672 + s26) * s26 - s12 * (s13 - t99 - t672) + s
     &16 * (s16 + ddreal(18) * s25) - t1117) * s26 + (s12 * (s16 + t110 
     &- t117) - t107 * (-s26 + t523 - t1118)) * s34 - (-(-t1121 * t26 + 
     &s16) * s45 + t26 * s34 * (-t61 * t1104 + s34 - t102 - t98) + t529 
     &- (t840 + t908) * s25 + (s13 + t99 - t103 + t106) * s26 + t1120) *
     & s45 + ((t670 - t805 - t114) * s25 - t116 * t282 + t16 * (t884 + t
     &316) + s12 * (t562 - t145 - t291 + t161 - t138) - s13 * (s16 + t29
     &3) + s16 * t890 + t1122 - t1123) * s46 - ((t145 + t1124 + t672) * 
     &s26 - (t110 + t113 + t1125 + t138) * s34 + (t16 * t1104 - t1127 + 
     &t292) * s45 + t123 + t1128 + t1129 + t1131 - t1130 + t136 + t528) 
     &* s56
      t1088 = ddreal(11) * t123
      t1089 = s25 - s16 - s26
      t1110 = (-s13 + t490 + s25) * s25
      t1121 = s25 + s16 + s26 - s34 - s45 + t113
      t1135 = t491 + t105 + t289
      t1136 = t116 * (-t121 + t258)
      t1137 = (t110 + t1105) * s26
      t285 = -s16 - t285 - t101
      t1138 = t26 * t98
      t1139 = t121 - t134
      t1140 = (t113 + t872) * s25
      t1141 = t1139 * t26 - (s12 + t104) * s16 + (s12 + t111 + t872) * s
     &26 + (s26 + t1138 - t105 + t109) * s34 + t1140
      t1142 = -t607 - t293
      t298 = -t26 * (t70 + t299) - t1053 - t1108 - (t488 + t841 + t1105)
     & * s26 + (-t117 + t826 + t697 + t101 + t908) * s34 + (s25 + t298 +
     & t238 + t872 - s45) * s45 + t162
      t1108 = -t300 - t490 + t144
      t1109 = -t26 * (-t1108 * s46 + t121 + t299) - t1060 - t1109 - (t48
     &8 + t826 + t908) * s26 + (-t117 + t840 + t1106 + t101) * s34 + (t6
     &01 - t293 + t872 + t289 - s45) * s45 - t828
      t1061 = (t1056 * s16 + s16 * t1061 + t70) * s25
      t1132 = (s13 * t140 + (t110 + t908) * s16 + t70 + t106 * t1132 + t
     &1111) * s26
      t1114 = (-t1114 - t102) * t123
      t1143 = (-t1089 - t293) * s45 - t26 * (t121 + t1086) - t394 - t830
     & - (s12 - t111 + t908) * s26 + (-s25 + t840 + t138 + t1105) * s34 
     &- t821
      t798 = t26 * t798 * s34
      t1144 = t798 - (s13 + t841 + t1134) * s26 - t761 - ddreal(10) * t1
     &72 - t684 - t361
      t479 = epinv * (-t1083 * s16 - (-s12 * t1135 - s34 * t285 - t1065 
     &- t1096 + t1136 - t1137) * s34 - (-s45 * t1142 - t1141) * s45 + s4
     &6 * t298 + s56 * t1109 + t148 * t1108 + t66 * t1108 - t1061 - t113
     &2 + t1114) + ((-s13 + t1056) * s16 + t26 * t828 - t506) * s16 + (t
     &1067 - t162 + t446 + t479) * s25 + ((ddreal(11) * s12 + s13) * s16
     & + (t283 + t102) * s25 - t506 + t70 + t1084 + t1088 + t1111) * s26
     & + t1144 * s34 + t1143 * s45 + (t26 * (-(s12 + t555 + t109) * s45 
     &+ t70) + t1053 + t830 + (s13 + t841 + t114 + t101) * s26 - (s25 + 
     &t360 + t826 + t159) * s34 + t1062 + t153 - t162 + t1110) * s46 + (
     &t26 * (s46 * t1121 + t121) + t1060 + t830 + (s13 + t826 + t159 + t
     &101) * s26 - (s25 + t118 + t840 + t159) * s34 + (s45 - s12 - t145 
     &+ t161 - t289) * s45 + t1049 + t1110 + t828) * s56 + t148 * t1121 
     &+ t66 * t1121
      t1062 = t26 * t135
      t1067 = s12 * (t360 + t159)
      t1083 = ddreal(12) * s16
      t1084 = t15 * t123
      t1110 = (t1083 + t118) * s26
      t1111 = t61 * s16 * (s16 + t98)
      t1121 = s16 + s34
      t1143 = t99 * s12
      t1144 = (-t117 + t106) * s16
      t1145 = -t134 + t162
      t1146 = t16 * t246
      t1147 = s16 - t117 + t1146
      t1148 = (t138 + t156) * s26
      t1149 = (t360 + t114) * s25
      t1150 = (-t872 - t289) * s13
      t1151 = -s16 + t117 - t256
      t1152 = t26 * t115
      t303 = -t246 + t303
      t1153 = s16 + s25 + s26 - s34
      t1154 = t61 * (t123 + t134)
      t1155 = (-t488 - t840 - t872 - s26) * s26 + (-t117 + t840 + t697 +
     & t114 + t101) * s34 + t1153 * s45 - t828 - t111 * t100 - t1143 - t
     &1154
      t292 = -t300 + t670 - t292
      t1156 = (t488 + t826 + t114) * s26
      t1157 = (-t117 + t840 + t697 + t101 + t908) * s34
      t1158 = (t601 - t293 + t99 + t289 - s45) * s45
      t1159 = t26 * (t121 + t299)
      t1160 = (-t1056 - t872) * s16
      t1161 = (t100 - t144) * s46
      t1162 = ((s25 + t145 + t117 + t106 + s26) * s26 + s16 * (t110 + t8
     &72) - t1117 + t70 + t106 * (s25 + t102)) * s26
      t1163 = -(t1121 * t26 + t246) * s45 - (s25 - t98) * s25 + (t488 + 
     &t113 + t109) * s26 + (-s25 + t145 + t118 + t840) * s34 - t394 + t1
     &111 + t497
      t1164 = -s25 * t247 - t1067 + t1080 - t1084 - t1110 - t361 + t778
      t1071 = t61 * (t123 + t121 + t134) + (t145 + t1071) * s26 - (s25 +
     & t906 + t840) * s34 - (s16 - t144) * s45 + t1107 + t828 + t1143 + 
     &t262
      t1054 = epinv * (-(-s12 * t1147 - s34 * t1151 - t1148 - t1149 - t1
     &150 - t125) * s34 - (-t303 * s45 - t26 * ((t98 - t1152 - s34) * s3
     &4 + t121) + t1111 - (s12 + t499) * s26 - t1140) * s45 - (-t1155 + 
     &t1161) * s46 - (-s46 * t292 - t1108 * s56 + t1060 + t1156 - t1157 
     &- t1158 + t1159 - t1160 + t828) * s56 - t1061 + t943 - t1162 + t10
     &34) + ((s12 + t99 + t101 + s26) * s26 - (-t145 - t111) * s25 + t10
     &93 - s12 * (s13 - t111 - t109) - t1133) * s26 + t1164 * s34 + t116
     &3 * s45 + t1071 * s46 + ((t26 * (s25 + s12 + s46) + t109 - t137 + 
     &t238 + s56) * s56 - t26 * t1145 + t1093 + t529 + (-s13 + t872 + t1
     &06 + s25) * s25 + (t110 - t117 + t159 + t106) * s26 + (s13 - t1091
     & - t1054) * s34 + (-s25 + s16 + s34 + s13 + t365 - s45) * s45 + (-
     &s13 + t1043 + t118 - t293 + t99 + s46) * s46 + t1144) * s56 + t111
     &3 - t943
      t609 = t26 * t609
      t1071 = t16 * t627
      t1091 = s16 + s13
      t825 = -s25 * t825
      t1107 = s16 * t1091
      t1111 = t26 * t153
      t1083 = t1147 * s34 - s12 * (t110 + t118 - t117) - s16 * (t102 - t
     &523) - (t1083 + t103 - t104) * s26 - ddreal(9) * t121 - ddreal(9) 
     &* t449
      t659 = (-(t98)**2 - t1099) * s16 + ((s13 + t109 + t841) * s16 + t3
     &61) * s25 + ((s25 + t117 + t99 + t106 + s26) * s26 + (s13 + t99 + 
     &t826) * s25 - t506 + t70 + t829 + t1093) * s26 + t1083 * s34 + (-(
     &s16 - t609) * s45 - (-s16 + t118) * s16 - t15 * s25 * t59 + t26 * 
     &(-s12 * (t281 + t111) + (s26 + s34 - s13) * s34 + t258) - t61 * s1
     &6 * t203 - t659) * s45 + (-(t246 + t144) * s46 + (-s26 - t57) * s2
     &6 - (t1043 - t117 + t159 + t697) * s34 - (-s12 + s16 + s25 + s26 -
     & s34) * s45 - t162 - t825 + t1107 + t1060) * s46 + ((-t26 * (s26 +
     & s34 + s46) + s12 + s13 - s16 - s25 - s45 - s56) * s56 - (-t1091 *
     & t61 + s25 - t840) * s26 - (-s13 + t840 + t697 + t114 + t101) * s3
     &4 - (t488 + t145 + t118 + t113 - t144) * s45 - (-s12 + s16 + s45 -
     & s13 + t521 + t111 + t289 + s46) * s46 + t70 + t822 - t825 + t1154
     & + t1111) * s56
      t825 = ddreal(17) * s16
      t829 = ddreal(14) * s12
      t1043 = t116 * t134
      t1083 = t2 * t121
      t1091 = ddreal(14) * s25
      t1099 = t26 * t456
      t1113 = t696 * s16
      t1133 = -s16 - t110 + t117 - t238
      t1143 = (-t116 * t98 - t109) * s16 + (s13 + t826 + t103 + t156 + t
     &501) * s26 - t26 * ((-t26 * t297 + s34 + t291 - t98) * s34 - t1068
     &)
      t1154 = -t293 - t1146 + t109
      t1163 = (s13 + t829) * s16
      t1164 = (t102 + t826) * s25
      t1165 = t26 * (-t153 + t162)
      t1106 = t16 * (t121 + t70) + t1043 + t1088 + (t829 + t825 + t104) 
     &* s26 - (t291 - t117 + t892 + t677 + t1113) * s34 - (-t521 + t840 
     &+ t1106 + t101) * s45 + t1163 + t1164 - t1165
      t1166 = -t1104 + t54 + t521 - t1087
      t1167 = s16 - s46
      t1168 = t129 * s34
      t1169 = t16 * t153
      t245 = -t1167 * t26 + t16 * t203 + t142 + t245
      t489 = ((s13 + t1064 + t1113 + t106) * s12 - t26 * (-s16 * (t110 +
     & t159) + t489)) * s26
      t1090 = (s13 * t1090 - (s12 - t523) * s16) * s16
      t1170 = t26 * (t121 * (s26 + s25 + t840 + t159 + t104) + t1076 + t
     &124)
      t1171 = (-s13 - t113 - t872 - t101 - s26) * t121
      t576 = (-t562 + t576 - t1099 + t840 + t109 + s56) * s56 - t1145 * 
     &t61 - t15 * t704 + t1082 + (-t523 + t841 + t908) * s25 + (-s13 + t
     &1091 + t897 + t156) * s26 + (s13 - t291 - t840 - t831) * s34 + (-s
     &25 + s12 - s26 + s13 - t99 + t144) * s45 + (t2 * t621 - s13 + s46 
     &- t142 + t145 - t521) * s46 + t70 + t114 * s12
      t621 = t798 - (t892 + t1105) * s16 - t1082 - t136 - t361 - (s13 + 
     &t111 + t841 + t825) * s26 + t549
      t704 = -t26 * t1115 * s16 + (-t117 - t840) * s13 + t70
      t798 = -(s16 + t293 + t105) * s45 - t16 * t971 - t26 * s16 * (s16 
     &- t1119) + (t16 * t98 + t159) * s25 + (t98 + t110 - t872) * s26 + 
     &(-t1085 * t16 + t360 + t908) * s34
      t674 = (t1153 * t61 - t54 + t840) * s46 + t1048 * t16 - t26 * (-t1
     &53 - t1049 + t162) + ddreal(10) * t123 + (t910 - t523) * s25 + (s1
     &3 + t829 + t1064 + t825) * s26 - (t101 + ddreal(13) * t107 + t892)
     & * s34 + (s25 - t840 - t674 - t1105) * s45 + t1043 + t1083
      t496 = epinv * ((t1133 * s34 + s16 * (ddreal(16) * s25 + t109) + (
     &t1118 + t697 + t825) * s26 + s12 * (s16 - t117 + t138 + t1064) - s
     &13 * (s26 + t872)) * s34 + (s45 * t1154 + t1143) * s45 + (s46 * t1
     &166 - t1106) * s46 + ((t245 + s56) * s56 + (t1044 - t727) * s45 + 
     &(-t145 + t588) * s46 + t116 * (t147 - t1168) - t26 * (s25 * t1167 
     &+ t121) + t61 * (-t70 + t316) + t696 * (s34 * t107 - t150) + ddrea
     &l(11) * t208 - ddreal(11) * t123 + s12 * (-t116 * t592 + t2 * t456
     & + s13 - ddreal(14) * s16 - t631) - s13 * (t565 + t118 + t109) + t
     &148 - t1169) * s56 - t489 + t1090 - t1170) + t704 * s16 + ((-s13 *
     & t15 + t106) * s16 - t840 * s13) * s25 + (t15 * t496 + t2 * t287 -
     & t26 * t506 + ddreal(12) * t123 + s12 * (s12 - s13 + t105 + t156))
     & * s26 + t621 * s34 + t798 * s45 + t674 * s46 + t576 * s56
      t506 = t1126 * s34
      t285 = (-s12 * t1135 - s34 * t285 - t1065 - t1096 + t1136 - t1137)
     & * s34 + (-s45 * t1142 - t1141) * s45 - s46 * t298 - s56 * t1109 -
     & t148 * t1108 - t66 * t1108 + t1061 + t1112 + t1132 - t1114
      t252 = (-s12 * t1094 - s34 * t252 + t1095 - t1096 - t1097) * s34 +
     & (-s45 * t839 - (t618 + t289) * s16 + t26 * (-s16 * t1126 - t121) 
     &+ t624 + s12 * (-t281 + t102 - t105) + t123 + t506) * s45 + t1063 
     &+ t1069 + t1076 + t1100 + t1101 + t1102 + t331
      t281 = s16 * (s16 + t1091)
      t298 = (t840 + t872) * s25
      t331 = (s12 + t102 + t618 + s26) * s26
      t565 = ddreal(14) * s26
      t576 = (t145 + t826) * s25
      t621 = ddreal(11) * s34
      t624 = (-s16 + t111 + t161) * s45 - t26 * ((t473 + t1138 - t105) *
     & s34 + t1140) + t394 - (s16 - s12 - t160) * s16 + (t98 - t110 - t8
     &72) * s26
      t156 = (s16 + t110 - t523 + t60) * s34 + s12 * (-t473 + t523 - t10
     &3) + (t138 + t908) * s13 + (-t565 - t1134) * s25 + (-t360 - t156) 
     &* s26 - t1093
      t473 = -t117 + t101
      t674 = s12 + s16 + s26 - s34 - s45 - t72
      t704 = t26 * t13
      t798 = ddreal(12) * t320
      t839 = t1 * t49 - t1008
      t910 = t1 * t997
      t971 = t274 * t165
      t528 = t971 * (t1 * t991 - t97) * (t26 * (t67 + t645) + (-s16 * t1
     &115 - t1040) * s16 + (t159 * s12 + s16 * t881 + t70) * s25 + (s12 
     &* (-s13 + t110 + t99) - t1117 + t281 + t331) * s26 + (t473 * s34 +
     & s12 * (-s16 + t117 - t105) + t107 * (-s26 + t523 - t672)) * s34 +
     & ((t115 - t555) * s45 + t26 * (-(s12 + t109) * s25 + (s34 - t98 - 
     &t289 + t1152) * s34) + t529 + (s13 - t291 + t99 + t106) * s26 + t1
     &120) * s45 + (-(-s16 - t1057) * s16 + (t872 + t841 + t103) * s26 -
     & (t291 - t117 + t1087 + t138) * s34 - (t100 + t101) * s45 + t1068 
     &- t162 + t1123 + t684) * s46 + ((t145 + t1124 + t103) * s26 - (t29
     &1 + t113 + t1125 + t138) * s34 + (-t61 * t820 - t1127 + t238) * s4
     &5 + t1068 + t1129 + t123 + t1131 - t1130 + t136 + t528) * s56)
      t68 = t971 * t68 * ((s12 * (-s13 + t102 + t103) - t1117 + t281 + t
     &331) * s26 + (-(s13 - t101) * s34 + s12 * (s13 - t291) + t107 * (-
     &s26 + t117 - t672)) * s34 + (-t169 + t26 * (t120 + t121) + (s13 + 
     &t490 - t618) * s26 + (-t98 + t291 - t102 - t289 + s34) * s34 + t12
     &3 - t105 * t167) * s45 + ((t110 + t840 + t109) * s26 - (-s13 + t49
     &0 + t1146 - s34) * s34 - (s16 + t101) * s45 + t824 + t828 + t298 +
     & t492 + t666) * s46 + ((t26 * (s12 + s16 - s34 + s46) + t238 + t58
     & + s56) * s56 + (-t1115 + s16) * s16 + (t110 - t117 + t99 + t106) 
     &* s26 + (-s12 - t491 - t1146 + s34) * s34 + (s16 - s34 + s13 - t10
     &1 + t365 - s45) * s45 + (-t61 * t203 - s13 + s46 + t118 + t499) * 
     &s46 + t298 - t1039 + t529) * s56 + t1033 - t258 * (s16 + t1092))
      t98 = t274 * s34
      t167 = -t1 * t463 + t52
      t59 = -t554 * (-t26 * (t67 + t1037 + t255) + t1034 + epinv * t895 
     &+ (s12 * t767 - t500) * s25 + (t788 + t932) * s26 + ((t114 + t289)
     & * s26 + t26 * t790 + (s12 + t99 + t289) * s12 + t449 - t796 + t79
     &7 + t1035) * s34 + ((s34 + t478 + t109) * s45 + t26 * (t821 + t134
     &) - t529 - (-t26 * t820 + s26) * s34 - t822 - t823 - t830) * s45 +
     & ((-t618 + t780) * s26 - t61 * t627 + s12 * (t562 + t521 - t831 - 
     &t101) + s16 * (s45 - t291 - t697 + t780) - t832 - t827) * s46 + (t
     &16 * (-s34 * t537 - t844) + t26 * (-t287 - t148 + t843 + t905) - t
     &61 * s16 * t59 + s12 * (s34 - t238 - t835 - t833) - t357 - t860 - 
     &t870 - t159 * s46 - ddreal(11) * t335 - t1036) * s56 + t943) + t63
     &3 * t119 + t700 * (t1054 - t1034 - t1062) + t926 * (t124 * t16 - t
     &1062 + t479)
      t119 = t362 * (s16 * t602 + (t601 + t102 + s26) * s26 - t170 + t20
     &9 - t390 + t471 + t604) - t386 * t629 + t422 * (t26 * ((-s16 - t64
     &7) * s25 - s26 * t83 + s34 * t573 - s56 * t573) + epinv * (t26 * t
     &564 + (-t548 + t555 + s16) * s16 + s25 * t676 + s26 * t435 + s34 *
     & t562 + t121 - t351 - t86 * t61) + s12 * (epinv * (t320 - t293) + 
     &t235 - t345)) + (t421 - t330) * (-s16 * t459 + (-t600 - t60) * s34
     & - (s16 + t144) * s45 + (t469 - t144) * s46 - t524 * t61)
      t169 = t8 * t1014
      t50 = t279 * (t1012 * (t10 * t119 + t21 * (t12 * t64 + t165 * t59 
     &+ t278 * (s25 * t280 + epinv * (s46 * t482 + (t514 + s56) * s56 - 
     &t152 - t290 + t518 - t519 + t142 * t246) - t288 * s46 + (t455 + t1
     &05 + s56) * s56 + t277 - t337) + t307 * (t61 * t515 * s25 + t250 +
     & epinv * (t329 + t323 - t296 + t305 + t317 - t290) + t309 - t337) 
     &+ t318 * (-t286 + t544 + t608 - t610 - t615 - t617) - t452 * (epin
     &v * (t329 - t544 - t608 + t610 + t615 + t617) + t26 * (t583 + t263
     & + t261 + t262)) - t484 * ((s16 - t256) * s45 - s46 * t482 - (t514
     & + s56) * s56 + t290 - t518 + t519) + t485 * (t286 + t323 - t296 +
     & t305 + t317 - t290)) + t8 * (t12 * t167 + t267 * t707)) + t699 * 
     &(-t190 * t21 * t1013 + t169 * (-t1 * t460 + t50)) + t465 * t10 * t
     &12)
      t50 = t50 + t1031 * (t10 * (t1012 * (t1015 * (t356 * t629 - t353 *
     & (epinv * t707 - t695)) + t1010 - t264) + t69 * (t72 - t704) * t10
     &15 + t971 * (t1 * t995 - t551) * t252 + t528 + t971 * (t910 - t768
     &) * t285 + t68 + t98 * t251 * t1014 * t699) + t21 * (t1 * (t978 * 
     &t699 * t1013 * t1015 + t1012 * t461) - t1012 * t51 + t971 * (-t559
     & * ((-s12 * t1147 - s34 * t1151 - t1148 - t1149 - t1150 - t125) * 
     &s34 + (-t303 * s45 + s12 * (-s26 - t1045 + t109) + t16 * s34 * t11
     &5 - t26 * (s26 * t247 - t506) + t61 * t1107 - t291 * s16) * s45 + 
     &(-t1155 + t1161) * s46 + (-s46 * t292 - t1108 * s56 + t1060 + t115
     &6 - t1157 - t1158 + t1159 - t1160 + t828) * s56 + t1061 + t1162 - 
     &t943 - t1034) - t571 * t285 - t574 * t252 - t622 * t895)) + t274 *
     & t18 * s34 * t699 * t1013 * t839)
      t59 = -t1024 + t669
      t64 = t1 * t977
      t68 = t1 * t989
      t115 = t1 * t1003
      t35 = t35 * t18
      t119 = -t275 * t558 - t447 * t508 + t59 * (t26 * (t124 + t1061 - t
     &1171) + (-t826 * s13 + (-t523 + t826) * s16 - t1038 * t26) * s16 +
     & (-t26 * (-s16 * (t872 + t103) + t1117) + s12 * (s12 - s13 + t1044
     & + t103)) * s26 + t156 * s34 + t624 * s45 + ((-t26 * t890 + s25 + 
     &t289 - t521 + t840) * s46 + t16 * t70 + t1042 + t1043 + t1083 - t1
     &165 + (s13 + t110 + t829 + t1113) * s26 - (-t523 + t892 + t1044 + 
     &t618 + t677) * s34 - (t495 * t61 + t138 + t840 + t908) * s45 + t11
     &63 + t576) * s46 + ((s25 - s45 - s13 + t345 + t118 + t1087 - t670 
     &+ s56) * s56 + (-s13 + t1116) * s16 + (-s13 + t110 + t1044 + t897)
     & * s26 - (t1056 - t160 + t1044 + t618 + t907) * s34 + (s12 - s26 +
     & s34 + s13 - t872 - t101) * s45 + (t2 * t346 - s13 + s46 + t111 - 
     &t142 + t145 - t621) * s46 + t576 + t70 - t896 + t1042 + t1043 + t1
     &36) * s56)
      t119 = t1012 * (t10 * t119 + t21 * (-t1079 * t616 - t575 * t659 + 
     &t799 * (-(t1133 * s34 - t26 * (-t1128 + t162) + t1093 + t818 + (s1
     &2 - t160) * s16 - (s13 - t1118 - t841 - t825) * s26) * s34 - (s45 
     &* t1154 + t1143) * s45 - (s46 * t1166 - t1106) * s46 - ((t245 + s5
     &6) * s56 - t61 * t70 - t1043 - t1088 - t1169 - t497 - (t829 + t104
     &) * s16 - (ddreal(12) * s12 + t160 + t1113) * s26 - (-t107 * t696 
     &+ s13 - t291 - t897) * s34 + (s13 + t826 + t1044 - t727 + t101 + t
     &157) * s45 + (-s13 - t145 + t111 + t588 + t772 - t897 + s46) * s46
     & - t1164 + t162) * s56 - t1090 + t489 + t1170))) * t165
      t91 = t91 * t21
      t125 = s34 * t462
      t3 = t1031 * (t21 * (t1015 * (t1012 * (-t244 * (epinv * t630 + t26
     &1 + t262 + t263 + t583) + t1016 * t762 * t1013) - t755 * t78) + t9
     &4) + t1012 * (t1015 * (-t12 * t986 - t240 * (epinv * t599 - t61 * 
     &t909 + s12 + s16 - s34 - t310 - t60) - t395 * (-t525 * epinv - t65
     & + t653 + t656 + t666 + t668 + t691) + t1 * t988 * t470 * t1014) -
     & t125) * t10 + t381 * (-t12 * t467 + t674 * t79 - t22 * (s25 + s12
     & + s16 - s34 + s56 - s13 + s46 + t60) * t1012) + t1015 * (t1013 * 
     &(t1027 * t21 * t1 + t1028 * t21 - epinv2 * (t3 * (t5 * ((propZ25 *
     & t25 * t15 + t2) * s25 + t302 + t798) - t19 * t819) + t20 * (t498 
     &+ t423)) * t1012 * t1014 - t91) - t1014 * t1032 * t10 - t1029 * t1
     &013 * t18 * s34) * t699)
      t5 = t443 * (-t26 * (t135 - t67 - t1037 - t255) - t1034 - epinv * 
     &t659 + (-s12 * t767 + t500) * s25 + (-t788 - t932) * s26 + (-t26 *
     & (t790 + t70) - t1084 - t1110 + t796 - t1067 - t449 + t1080) * s34
     & + (-(t609 + t109) * s45 + t1139 * t16 - t26 * t821 + t830 + (t16 
     &* t346 - t111 + t872) * s34 + t822 + t823) * s45 + ((-t588 + t618)
     & * s26 + t827 + t832 - s12 * (t562 + t670 - t831 - t101) + s16 * (
     &-s45 + t291 - t588 + t697) + t1071) * s46 + (s12 * (-t293 + t238 +
     & t835 + t833) + (t451 + t109) * s16 + (ddreal(11) * s25 - t521) * 
     &s26 - t16 * ((s25 + s16 + s46) * s34 - t844) - t26 * (-t134 - t287
     & - t148 + t843 + t905) + t860 + t870 + t1036) * s56 - t943) - t733
     & * (t1098 * t16 - t26 * (t135 + t1171) + t496)
      t19 = -t192 * t585 - t214 * t599 - t216 * (epinv * t585 - t584 * t
     &61 + s12 - s34 - s46 - t310) - t229 * t9 + t32 * (-epinv * t229 - 
     &s13 + s25 + s46 + s56 + t310 - t596 - t60) + t342 * t525 - t347 * 
     &((-t367 - s56) * s56 + t172 + t173 - t516 - t425 + t394) + t404 * 
     &(epinv * (-(-t367 - s56) * s56 - t172 - t173 + t516 + t425 - t394)
     & + (-t16 * t58 - t60) * s34 + (-t297 - t161) * s45 + (t340 - t161 
     &+ t365) * s46 + (t42 - t161 + s56) * s56 + t434 + t440 + t446)
      t5 = t321 * (t10 * t19 + t18 * (-t222 * (epinv * t586 + t26 * (-s1
     &3 - s16 + s25 + s45 + s56)) - t33 * (epinv * t229 + t695)) + t21 *
     & (-t14 * t229 + t165 * t5 - t197 * t586 - t215 * t599 - t236 * (ep
     &inv * (-t838 + t911 + t142) - t239 - t698) - t370 * (-s46 * t532 -
     & s56 * t511 + t123 + t535 - t540 - t542 - t546) + t429 * ((t105 + 
     &t289) * s25 + epinv * (s46 * t532 + s56 * t511 - t123 - t535 + t54
     &0 + t542 + t546) + (t455 + t291 + s56) * s56 - t549 + t550 + t552)
     & - t81 * t630 - t742 * ((-t112 + s26) * s26 + (s56 + s45 - s13 + s
     &46 + t478) * s56 + t246 * t84 + t644 + t761) * t1013) + (-t30 - t1
     &78 - t180) * t55 * t8)
      t19 = t807 * t760
      t20 = t972 * ((t100 - t756) * s25 + t171 - t626 + t758)
      t42 = t561 * (t26 * (-t124 + t135 - t1037) + epinv * t1079 + ((-s1
     &6 + t104) * s16 - t113 * t476) * s25 - (t16 * s25 * t132 - (s13 - 
     &s16 + t101) * s13 + t70 + t1042 + t113 * (s25 + t99) + t1072) * s2
     &6 - (-(t159 + t897) * s16 - (t488 + t1044 + t841 + t118) * s26 - t
     &26 * (t761 + t70) + t1080) * s34 - (-(-s26 + t1045) * s45 + (t727 
     &+ t101) * s26 - t16 * t472 - t26 * ((-s12 + s34 + s25) * s25 + t12
     &1 + t123 - t778) - s12 * s26 - s16 * (t360 - t780 - t101)) * s45 -
     & (t1048 * t26 + t797 + (s13 + t110 + t1046) * s26 - (t906 + t111 +
     & t826) * s34 - t162 + t1049 + t1050 - t1051 + t1053) * s46 - (t26 
     &* (s46 * t1055 + t121 + t287) + t797 + (t132 * t15 + s13 + t826) *
     & s26 - (t145 + t118 + t111 + t840) * s34 + t1050 + t828 + t1058 - 
     &t1059 + t1060) * s56 - t148 * t1055 - t66 * t1055 - t1081)
      t38 = t38 * s34
      t42 = t279 * (t10 * (t1014 * (t1012 * (-t166 * ((s34 - t168 + t756
     &) * s34 - t170 + t171 + t172 + t173) - t468 * t470 - t699 * t96) +
     & t339 * t699) - t1012 * (t660 * (epinv * t275 - (t1066 * s34 + t26
     & * (-(t872 + t106) * s25 + t162 - t836) - t61 * ((t283 + t837) * s
     &26 + t123) - t684) * s34 - (t457 + t1073) * s45 - t1070 - t1074 - 
     &t1075 - t1076 - t1077 - t1078 - t1081 - t783) + t845 * t1047) * t1
     &65 - t72 * t1007 * t1012) + t8 * t1012 * (-t1014 * t699 * t1021 + 
     &(-t31 - t992 - t993) * t55 * epinv - t38 * t755) + t1012 * (t1013 
     &* (t699 * t982 - t762 * t975 + t19 + t20) + t374 - t375 - t42 * t1
     &65) * t21)
      t55 = -t358 + t56
      t56 = t1012 * t10
      t72 = t1 * intHs16s25s26s34s56x1411D6eps1()
      t96 = intHLs16s25s26s34s56x1114D6eps0() * t21
      t65 = t1031 * (t1015 * (t21 * (-t1012 * (t269 * (epinv * (s34 * t5
     &77 + s45 * t48 - t204 - t579 - t580) + t261 + t262 + t263 + t583) 
     &+ t483 * (s45 * t607 - (s12 + s16 - s34 - s45) * s56 + t259 - t299
     & - t608 - t273)) - t466 * t755) + t1012 * (t12 * (t10 * t328 - t10
     &05 * t21) + t383 * t674) * s34 + t56 * t55 * (t668 + t691 - t65 + 
     &t656 + t666 + t653)) + t8 * (t72 - t270) + t96)
      t104 = ddreal(1) / ddreal(18)
      t3 = t163 * t1031 * (t21 * (-t73 * t755 * t1015 + t95) + t98 * (-t
     &12 * t21 * t1006 + t18 * t984 * t674)) + t196 * t50 - t104 * (t42 
     &+ t5 + t3 + t279 * (t1012 * (-t115 * t21 * t760 * t1013 + t1014 * 
     &(t68 - t475) * ((-s25 - s16 + s34 - s45 + s13 - s46 - t60 - s56) *
     & s56 + t477 * t84 + t258 - t299 - t304 - t480) * t10 + t169 * t174
     & * (t177 * t699 + t176) + t35 * s34 * t674) + t699 * (t1013 * (t10
     &12 * (t21 * (-t979 + t980 + t981) + t186 * t18) + t64 * t21) + t10
     &12 * (-t257 * t8 + (-t1 * t37 + t28 + t29) * t8 * s34 + t276 * t10
     &) * t1014) + t119)) - t65 / ddreal(3)
      t5 = t1 * t11
      t12 = t578 * t976 - t368
      t28 = t11 * t82
      t29 = s12 * (s56 + t388 - t144) + t28
      t37 = -t26 * t459
      t42 = t924 - t37
      t50 = s12 + s16 + s26 + s46
      t65 = t82 + t113
      t98 = -s13 - s16 + s25 + s34 - s46 + s56
      t110 = s16 + s26 - s34 - s45 - s46
      t119 = s13 * t110
      t132 = t248 * s34
      t136 = t779 * s45
      t145 = t26 * (t121 + t148 + t338)
      t156 = t98 * s12
      t159 = (t116 * t914 - t1093 + s16 * (s25 - t427 - t118 + t161) + t
     &335 - t119 - t132 + t136 - t145 + t156) * s12 - t510 * ((s16 - s34
     & + t312) * s16 + t143 * t348)
      t160 = t116 * t673
      t166 = -t678
      t169 = t86 * t15
      t170 = t16 * t654
      t172 = -s13 + s25 - s46 + s56
      t173 = t441 * s34
      t176 = (-s12 * t172 + (t118 - t144) * s46 - t26 * (s26 * t531 - t1
     &48 - t338) + s16 * (t548 + t212 + t60) + t119 - t335 + t951 + t173
     &) * s12 + t11 * t143 * t510
      t110 = t751 * t110
      t229 = t747 * s26
      t98 = t36 * t98
      t245 = t242 - t417 - t405
      t250 = t84 * s25
      t251 = (s25 + s34 + s45) * s34
      t252 = (-s46 + t950) * s45
      t255 = s45 - s46 - t405
      t256 = s46 + t406
      t257 = (-t242 - t417) * s56
      t263 = s25 + s45 - s56
      t264 = t11 * t509 * t510
      t193 = (s12 * t193 - t26 * (t263 * s34 + t121 + t134 - t335 + t966
     &) - t61 * (t366 + t123) + s16 * (s45 - t160 - t412 + t111) - t119 
     &+ t151 - t324 + t86 * t116) * s12 + t264
      t4 = t41 * (t61 * (t834 * t184 + t7) + t26) + t7 * (t4 + t6 * (t83
     &4 * t26 + ddreal(1)))
      t6 = -s13 + s25 + s26 - s46 + s56
      t7 = s34 - t414
      t41 = t144 * t143
      t184 = s16 * (t976 - t393 - t142)
      t267 = t393 * s34
      t275 = s25 * t578
      t276 = t61 * t714
      t280 = t237 * t320
      t281 = t426 * s46
      t285 = s34 - t177
      t286 = (t116 * s46 * (t405 + s26) + t16 * (t325 - t567) - t26 * (t
     &281 + t280 - t720) - t61 * (-t228 * t456 - t173 - t316) + (-s34 - 
     &t418 - t848) * s26 + s34 * (-epinv * t441 + s34) + s45 * t285 + t1
     & * (t119 - t184) - t736 - t36 * t6) * s12 + t510 * (-t220 * t7 + t
     &26 * (t325 + t879 + t275 + t524) - t877 + t276)
      t287 = -s13 + s25 + s26 + s34 + s56
      t288 = -t84 + t293
      t290 = s16 * (t48 - t414 + t161)
      t292 = epinv * t320
      t296 = t320 * s46
      t89 = t89 * epinv
      t297 = s34 * t456
      t276 = (t116 * t714 - t26 * ((-s25 + t292) * s45 - t231 - t296) + 
     &t61 * (t560 + t297) - t885 - s25 * t407 - (-s56 + t855) * s34 - (-
     &s25 - s45 + t89) * s46 - t1 * (-t119 + t290) - t415 - t715 + t720 
     &- t736 - t36 * t287) * s12 + t510 * (-t220 * t288 + t275 + t276 + 
     &t325 + t524 + t879 - t880)
      t298 = t26 * (s25 - s13 - s46) + s26 + t412
      t302 = -s26 + s34 + s46
      t303 = (-s12 * t298 + (t212 + t805 - t780) * s46 - t116 * t541 + t
     &26 * (-t335 + t121 + t951 + t119 + t231) + t61 * t302 * s56 + (t34
     &3 * t61 + s16 - s34 - t111 + t158 - t221) * s16 + t208 - t66) * s1
     &2 + t510 * (t312 * t348 + s16 * (s16 - t127 + t212 + t289))
      t304 = s16 * (-t26 * t811 + s25 + t118 + t144 + t221)
      t305 = t61 * (-t325 + t121 + t654)
      t309 = t26 * (t123 + t891)
      t317 = t406 + t405
      t323 = t177 * s26
      t329 = s56 + s46 + t405
      t331 = t264 * epinv
      t263 = (-(-s26 - s56 - t398 - t220 + t418 + t241 - t751) * s12 + t
     &116 * t329 * s26 - t16 * t813 - t26 * ((epinv * t263 + s56) * s34 
     &- t148 + t415 - t66 + t852 - t323 + t237 * t311) - t61 * ((-s26 + 
     &t397) * s26 + t151) - (t854 * s16 - t116 * t317 - t61 * t940 - t89
     &8 - t925 + t939) * s16 - (-s25 - s45 - t847) * s34 + t110 + t134 -
     & t230 - t769 - t862) * s12 + t331
      t311 = -t141 + t75
      t172 = (-t36 * t172 - t26 * ((s25 - t405) * s26 + t257 - t720 - t7
     &36) + t61 * (t322 - t913 - t86) - t722 + (t16 * t417 + t26 * t735 
     &+ s16 - t237 + t289 + t413 - t848) * s16 + (s45 + t848 - t177) * s
     &34 + (s25 - t242) * s46 + t110 + t134 + t933) * s12 + t510 * ((-s2
     &5 - s45 - t956 + t417 - t293 + s16) * s16 + s25 * t311 - t379 * (-
     &t1 * t143 + s34))
      t289 = -t116 + t408
      t337 = t592 * s26
      t285 = (-t36 * t298 + s46 * (t405 * t2 - t960) - t116 * (-s56 * t9
     &49 - t402) + t15 * t322 - t16 * (-t513 - t720 - t620) + t26 * (s34
     & * t285 + s45 * t928 + t110 + t323 - t852) + t61 * (-t397 * t302 +
     & t209 - t337) - ddreal(7) * t556 - (-(ddreal(9) - t730) * s46 - t2
     &89 * s56 - t26 * (s26 + s45 + t177) + t61 * (s25 + s34 + t752) - t
     &405 - t220) * s16 + t121 - t721 - t919) * s12 - t510 * (-t16 * t71
     &4 - t26 * ((-s26 * t857 + t147) * epinv + t134 + t454 + t545) + t9
     &53 - (-t212 * epinv - s25 - s56 + t220 + t60 + t734 - t737) * s16 
     &+ s25 * (s26 - t62 - t432) - t121 + t582)
      t282 = t26 * t510 * (s16 * (s45 + s46 - t241) + s34 * (s34 + s56 -
     & t406) + t275 + t560 + t663) - s12 * (s12 * (-epinv * (s25 + s16 +
     & s34 + s46 + t423) + s25 + s26 + s34 - s46 + s56 - t751) - t116 * 
     &t481 + t16 * (s26 * t127 - t338 + t884) - t26 * (t664 + t524) - t6
     &1 * ((s25 + s16 + s34) * s45 + t132) - epinv * (-t917 + t304 + t55
     &3 - t627 + t955 + t309 + t305) + s56 * t107 - t110 + t282 + t668)
      t298 = t26 * (t793 + t134 + t350) - t885 + s12 * t355 + t121 + t12
     &3 - t211 - t513 - t650
      t302 = s16 * t84
      t338 = t578 * s56
      t339 = s34 * t82 + t341
      t340 = s12 * t364
      t341 = t578 * t25
      t346 = s12 * t373 + t341
      t348 = (s16 - t924 + t60) * s16 + t143 * t379 - t382
      t350 = t369 * epinv
      t313 = -s26 * t11 - t313 + t314
      t353 = -t26 * s46 * t203 + s16 * (t343 + t345) + (-t811 - t293 + s
     &26) * s26 - t147 - t314
      t355 = (s12 + s16 + s26 - s46 - t75) * s56
      t198 = t198 + t345
      t358 = s12 + s16 + s26 + s34 + s45 + t345
      t364 = epinv * t358
      t365 = s16 + s26 - s46
      t367 = t365 - t75
      t131 = (s26 - t140) * s26 + (s25 - s26 - s13) * s34 + (s16 - s34 +
     & s45 - s13 + t60 + s46) * s46 + (t131 + t312 + s56) * s56 + t253 -
     & t299 + t507
      t57 = t26 * (t123 + t66 + t148) + t492 + (-t57 + t872) * s26 + (t5
     &7 - t238 - t109 + s34) * s34 + s45 * t308 + (t283 + t118 + t99 - t
     &144) * s46 + (t146 * t16 + s12 + t118 - t293) * s56 + t120 + t162 
     &- t299
      t120 = epinv * t11
      t140 = t342 * t346
      t146 = epinv * t79
      t186 = t186 * t699 * t1013
      t253 = t742 * (t48 * t578 + t644)
      t160 = -t508 * t193 - t558 * ((-t116 * t166 - t26 * (s45 * t652 + 
     &t121 + t132 + t148 + t513) - t1093 + s16 * (s25 - t427 - t160) - t
     &119 + t335 - t649 + t169 + t156) * s12 - t510 * (t26 * (s16 * t676
     & - t151) + (s26 - t811 - t144) * s26 + t123 + t134 + t271 + t316 +
     & t649)) + t660 * ((t116 * s46 * (-t405 + t229) + t26 * ((epinv * t
     &652 - s56) * s45 + t415 - t717 - t720 - t736 + t851) + t170 - t420
     & - s16 * (t116 * (s46 + t406 + t405) - t161 - t729 - t961) + (-s25
     & - s34 + t417) * s45 - t132 + t745 + t965 - t110 - t1093 * t1 + t9
     &8) * s12 + t510 * ((t649 + t134 + t271 + t316) * epinv + s16 * (-t
     &26 * (t406 + t708) + t924) + s26 * (s34 + t242 - t417 - t702) + t5
     &16 + t705 - t405 * t345)) - t845 * t263
      t193 = s12 * t319
      t28 = t1012 * (t18 * (t443 * t282 - t561 * ((-t16 * s34 * (t228 + 
     &t923) + t26 * ((-s25 - t848) * s46 - t716 - t720 - t736) - t61 * (
     &t251 + t410 + t316) + t169 + t722 + s16 * (-t116 * t743 + t242 + t
     &418 - t915) - t110 + t965 - t252 + t405 * t248 + t98) * s12 + t510
     & * (-t26 * (t251 + t250) - s16 * (-t245 + t591 - t144 + t220) + s2
     &6 * (t245 + t144) - (-t241 - t242) * s46 - t736)) + t575 * (t293 *
     & t82 * t510 + s12 * (s12 * (s25 + s16 + s34 - s13 + s46 + t423) + 
     &t66 - t119 + t304 + t553 - t694 - t627 + t955 + t309 + t305)) - t6
     &16 * ((-s16 * (t128 - t505 + t161) + t119 + t132 - t136 - t335 + t
     &557 + t145 - t170 + t1093 - t156) * s12 + t28 * t510) - t733 * t28
     &5 + t799 * t303) + t8 * t160) * t165
      t28 = t279 * (t1013 * (t21 * (t1012 * (t699 * (t27 * t999 - t979 +
     & t980 + t981) - t115 * t760) + t64 * t699) + t1012 * (t1017 * t12 
     &- t253) * t18 - epinv2 * t699 * t4 * t1012) - t193 * t81 * t21 * t
     &1012 + t56 * t165 * t59 * ((-(t26 * (s25 - s16 + s34 - s13 - s46) 
     &- t786) * s12 + s26 * (t360 - t588) + (t212 - t780 + t157) * s46 +
     & t116 * (t123 + t513) + t26 * (-t335 + t132 - t136 + t119) + t61 *
     & t757 - s16 * (t26 * t810 + s45 - t157 - t158 + t727) - t208 - t65
     &5) * s12 + t11 * (s16 + t312) * t510) + t28 - t1012 * (t10 * t1007
     & + t8 * t993) * t578 * epinv)
      t59 = t699 * t1019
      t115 = t8 * t467
      t19 = t1031 * (t1012 * (t10 * (t1015 * (-t11 * (t30 + t178) - t192
     & * (s12 - s16 - s26 + s34 + s45) - t216 * (s12 * t27 + t5) - t347 
     &* t298) - t125) + t1015 * (t21 * (s12 * t53 + t11 * t14) + t8 * (-
     &(t174 + t475) * s12 - t11 * t9)) + t1015 * (t21 * (t699 * (t1 * t1
     &018 + t1000 * t27 + t982) + t19 + t20) - t59 * t18 * s34) * t1013)
     & - t1015 * (t21 * t78 + t115) * t578 - t91 * t699 * t1015 * t1013)
      t19 = t19 + t28 + t279 * (t8 * (t578 * (-t1012 * t180 - t146) + t1
     &012 * (s12 * t989 + t11 * t32) * t1 - t22 * t50 * t1012) + t56 * (
     &-t214 * t358 - t240 * (-s34 - s45 + t364 - t62 - t37) - t395 * (s1
     &2 * (s16 + s26 + s56 + t398 - t241) + t341 * epinv) - t404 * (-t26
     & * (s26 * t100 + t653) + epinv * t298 + s34 * t842 + (-t11 - t113)
     & * s46 - t121 - t123 - t314 + t844) - t578 * t986 - t120 * (t31 + 
     &t992) + t38 * t40 - t140) + t1012 * (-t348 * t370 + t186) * t18 + 
     &t1012 * ((-s12 * t464 + t139) * s34 - t197 * t198 - t215 * t358 - 
     &t222 * (epinv * t198 - t698) - t236 * (t364 - t704) - t244 * (t268
     & * t319 + t231 + t261 + t262 + t338) + t429 * (-epinv * t348 + t23
     &1 * t16 - t355 + t550 + t552 - t668) + t33 * t34 * t11) * t21)
      t6 = -t176 * t574 + t554 * t286 + t559 * ((s12 * t287 + s45 * t218
     & + t116 * t654 - t119 + t121 - t132 + t290 - t326 - t594 - t625) *
     & s12 + t510 * (s16 * t288 + t41)) + t622 * ((s12 * t6 - t26 * (t32
     &5 + t148 + t513) + t61 * t888 - t119 + t121 - t173 + t184 - t316 +
     & t442 + t267) * s12 + t510 * (s16 * t7 + t41)) + t633 * t172 + t70
     &0 * t276 - t926 * (((t360 + t212) * s34 + t116 * (s46 * (-t406 - t
     &405) - t402) + t26 * (t257 - t322 - t720 - t736) - t61 * (t410 + t
     &627) - (-t134 - t231) * epinv + s16 * (-t116 * t256 - t15 * t918 -
     & t409 + t961) - t110 - t209 - t252 + t965 + t98) * s12 + t510 * (-
     &t26 * (s34 + s45 + t847) * s34 - (-t26 * (-t406 + t847) - t144 - t
     &255 + t220) * s16 + s26 * ((t61 - t710) * s34 + t242 - t417) - (s2
     &5 - s34 - t242) * s46 - t322 - t736))
      t7 = s34 * t984
      t20 = t987 * t578
      t28 = t63 * t578
      t34 = s34 * t1006
      t37 = (t249 - t190) * t21
      t38 = t1031 * (t1012 * (t1015 * (t18 * (t571 * t159 * t165 - t318 
     &* t339) + t8 * (s12 * t52 - t356 * t29)) + t21 * (t1015 * (-t452 *
     & (t26 * (t261 + t262 + t338 + t231) + epinv * t339) - t28 - t34 * 
     &t65) - t51 + t265) + t1026 * ((t42 * t983 + t902 * (epinv * t42 + 
     &s16 + s26 + s46 + t113)) * s34 + t165 * (-s12 * (t131 * t486 + t57
     & * t97) + t1 * (s12 * t991 * t57 - t159 * t997 + t176 * t995)))) +
     & t73 * t21 * t367 * t1015 + t381 * (epinv * t69 + t465) * t11 + t1
     &015 * (-t18 * s34 * t1008 * t1012 + t37) * t1013 * t699)
      t6 = t38 + t321 * (t10 * (t362 * (t369 + t340) + t422 * (s12 * (t2
     &6 * (-s56 - s46 + t241) + t424) + t350)) + t165 * (t10 * (t36 * t9
     &90 * t131 + t159 * t768 - t176 * t551) + t18 * t6) + t21 * (-t278 
     &* (epinv * t353 - t277 - t294 + t306 + t355 + t668) - t307 * (-epi
     &nv * t313 + t26 * t668 - t322 * t61 + t171 - t315) + t313 * t485 -
     & t484 * t353 - t7 * (epinv * t65 - t704) + t20) + t8 * (t29 * t386
     & + t428 * (s12 * (t82 - t144) + t332 * t82) - t36 * t463) + t241 *
     & t49 * t18 * t699 * t1013)
      t6 = -t126 * t403 * propW16 * s16 * t185 * t699 * t1012 * t88 + t1
     &96 * t6 - t104 * t19 - t403 * (t1023 * (t21 * (t1012 * ((-t1005 * 
     &t50 + t430 * (-epinv * t50 + s12 + s16 + s26 - s34 - s45)) * s34 -
     & t269 * (epinv * (s46 * t203 - t302 + t86) + t231 + t261 + t262 + 
     &t338) - t483 * (-t152 - t273 + t86)) + t11 * t466) + t385 * (s16 +
     & s26 - s34 + s46 + s56) * t8 * s12 + t56 * s34 * t50 * t328) + t10
     &12 * t219 * t699 * t185 * propW16) / ddreal(3)
      t19 = t11 * t578
      t29 = (s56 + t698 + s12) * s12 - t19
      t36 = s16 + s26 - s34 - s45 + s56
      t38 = -s12 * t36 + t19
      t40 = t578 * t319
      t41 = s16 + s26 - s56
      t42 = s12 * (t41 - t345 - t589) + t201 * t578 + t361
      t5 = t5 * t578
      t49 = s12 - s16 - s26 - s45
      t50 = s12 - s34 - s45 - s46 - s56
      t52 = epinv * t50
      t11 = s12 * (t11 * t61 + s56) - t19 + t361
      t57 = t26 * t38
      t63 = t373 * s56
      t65 = t520 + t522
      t78 = t26 * (s34 * t578 - t351) + (s16 + t65) * s16 + (s26 + t520 
     &- t75) * s26 - t63
      t82 = s25 + s26 - s34 + s56
      t91 = (s16 * t82 + t380 * t878 + t260 + t553) * s12 + t40 * s34
      t98 = (t380 * (t893 - t293) + t204 + t553 + t579) * s12 + t578 * (
     &t595 - t453)
      t110 = -t345 * s34 + t275 + t302 + t325 + t894
      t119 = s12 * s46
      t131 = t578 * (t391 + t325 + t879 + t453 + t87)
      t136 = t101 * t379
      t139 = t19 * t380
      t65 = ((t380 + t101) * s12 + s16 * (s46 + t423 + t101) + t380 * t6
     &5 + t136) * s12 - t139
      t145 = t384 * t578
      t152 = t111 * t70 + s12 * (t111 * t379 + (s16 - t924 + t628 + t111
     &) * s16 + t380 * (-t74 + t423)) - t145
      t156 = s34 - t397 - t417
      t159 = t61 * (t242 - t405)
      t160 = t405 * t74
      t169 = t26 * t319
      t136 = ((t512 + t101) * s12 + s16 * (t48 * t61 + s45 + t345) + t38
     &0 * (t300 * t61 + s56 - t161 - t345) + t136) * s12 - t578 * (s16 *
     & (t343 + t200) + t782)
      t170 = (t110 + t119) * s12
      t172 = s13 * (s16 + s26 + s34)
      t173 = t26 * t301
      t58 = -t26 * ((-t299 - t1086) * s16 + t149 + t67) + ((t494 + t106)
     & * s26 + t99 * t58 + s12 * (s12 + t494 + t872)) * s26 + (-s12 * t2
     &47 - t107 * t473 + t231) * s34 + (-(s16 + s26 - s34 - s13) * s45 -
     & t26 * (s16 * t248 + t121 + t517) - s12 * t137 - t123 + t132 + t17
     &2 - t150 * t61) * s45 + ((s12 + s34 + s13 - t54 - t118 - t99) * s4
     &6 - t26 * (t121 - t299 - t828) - (s13 - t826 + s16) * s16 + (-s13 
     &+ t1056 + t111 - t109) * s26 - (t107 + t494 + t106 - s34) * s34 + 
     &(s34 - t360 - t173 - t872) * s45) * s46 + ((s12 - t427 - t54 - t11
     &8 - t99) * s56 - t148 * t15 - t26 * (t123 - t299 + t162) - t492 + 
     &(t494 + t840 - t872) * s26 + (-t107 - t173 + s34) * s34 + (-s45 - 
     &t494 - t238 - t109) * s45 + (t137 + t113 - t310 - t697 - t114) * s
     &46 + t1144 + t70) * s56 + t653 * t487
      t106 = t26 * t343
      t99 = (-(s12 - t117 + t105) * s26 - s12 * (s16 + t291 - t117) + t9
     &9 * t112) * s26 - (-t101 * t168 + t117 * t107) * s34 + (-(s25 + s3
     &4 + s13 - t164) * s45 - t61 * t1089 * s34 + t116 * t668 + s12 * (t
     &1089 * t26 + t137) - s16 * t107 - t134 - t172) * s45 + ((t503 + t1
     &06 + s46) * s46 - (-t1041 + s16) * s16 - (t476 + t840 + t101) * s2
     &6 - (s34 - t26 * (s12 + s16 + s26 - s13) - t101) * s34 - (t117 - t
     &238 - t837 - s45) * s45 - t70 + t1039 - t101 * t100) * s46 + ((-s1
     &2 + s16 + s34 + t62 + t106 + s56) * s56 + t26 * (t153 + t162) + t6
     &1 * (t148 - t299) + (s16 - t284 - t101 + s26) * s26 + (t26 * (s16 
     &+ s26 - s13) + s12 + t101 - s34) * s34 - (t26 * (s12 - s26 + s13) 
     &- t101 - t1121) * s45 - (-t16 * t343 + t26 * t493 - t1126) * s46 -
     & t836) * s56 - t299 * (s12 + t99)
      t99 = t26 * (-t258 * t100 + t566) - t99
      t101 = s16 * (-s46 + t111)
      t106 = t379 * s46
      t109 = -s16 * t431
      t112 = t26 * t166
      t113 = s16 * (t592 - t345)
      t117 = -t26 * s46 * t141
      t118 = (t161 + t142) * s46
      t132 = t345 * t70
      t137 = t396 * s16
      t150 = -(s25 - s12 - s16 - s26 + s34 + s45) * s56 - t372 + t530 + 
     &t758
      t162 = t48 + t345
      t166 = s45 + t169
      t168 = t116 * t48
      t172 = s25 * t657
      t173 = t66 - t172
      t176 = s13 * t166
      t184 = t2 * t210
      t190 = t16 * t48
      t198 = t116 * t84
      t201 = ddreal(12) * s34
      t203 = s25 * (t243 - t670 - t772 + t501)
      t204 = t696 * s56
      t218 = ddreal(14) * s46
      t219 = t26 * t387
      t245 = t811 + t219
      t219 = t811 - t238 + t219
      t211 = t26 * (-t211 + t123 + t121 - t650)
      t247 = t209 * t379
      t248 = t207 * t379
      t249 = s25 - s45 + s46
      t251 = -s25 + s26 - s34 + s56
      t252 = s13 * (s16 * t219 + s26 * t245 - t211 + t594 + t789)
      t257 = ddreal(17) * s56
      t258 = t766 * t378
      t260 = t162 * s12
      t261 = t363 + t345
      t262 = t416 - t86
      t263 = -t116 * t435
      t268 = t61 * t435 + s46
      t273 = -t26 * t637 - t144 + t238
      t275 = s25 * t273
      t276 = ddreal(10) * t121 - ddreal(10) * t351
      t277 = t366 * s46
      t282 = t445 * t153
      t283 = t641 * s34
      t284 = (t335 + t864 + t134 + t133) * s26
      t285 = t207 * t74
      t286 = (t153 - t66 + t553 - t524 - t85) * s46
      t287 = t147 - t121
      t288 = ddreal(12) * s56
      t290 = (-t816 - t288) * s45
      t291 = t261 * s12
      t294 = t61 * t84
      t298 = t121 - t524
      t299 = t128 * t134
      t300 = t637 * s25
      t221 = ((-s12 * t363 + t16 * t325 + t26 * (t335 + t336 - t123 + t6
     &6 + t334) + t61 * t296 + s16 * (-t26 * t1104 + t310 + t373) + (-s2
     &5 + s26 - s45 - s46) * s34 + t148 - t176) * s12 - t116 * ((t147 - 
     &t366 + t86) * s56 - t444) - t26 * ((-t581 - t134 + t300 - t148 + t
     &351 + t590) * s26 + s45 * t945 + t614 + t690 - t754) - t61 * ((t33
     &5 + t66) * s34 + s45 * t863 - t643 + t748) + ((t221 + t478 + t294)
     & * s16 + (t359 + t451) * s56 + t16 * t298 - t26 * (-t357 - t148) -
     & t792 + s26 * (-s34 + t310 + t791 + t803) + t275 - t516) * s16 + (
     &s45 * t333 - t121 + t148) * s34 + t252 + t67 + t686 + t775 + t299 
     &+ t366 * (t360 + t451)) * s12 + t40 * t293 * t510
      t48 = t26 * t48
      t301 = s25 - s34 + s45
      t302 = s46 * t520
      t303 = t231 * t74
      t304 = t231 * t199
      t305 = s26 - s34
      t306 = s46 * t305 - t524
      t82 = ((t26 * (s26 * t301 + t121 + t209 + t448 + t66) - t61 * t694
     & + t557 - t771 + s16 * (-s34 + t212 + t48 + t142) + t148 - t176 - 
     &t271 + t119) * s12 - t116 * (t147 * t878 - t86 * t531) - t15 * s34
     & * t306 - t16 * s56 * ((-s45 + s56) * s34 - t121 + t208) + t26 * (
     &(t581 + t209 - t153 + t66 + t207 + t661) * s26 + s45 * t777 - t686
     & - t688 - t304) - t61 * (t948 + t685 + t536 - t444) + ((t26 * t82 
     &+ t294) * s16 + (-t451 - t658 + t144) * s34 + t16 * (t121 - t153) 
     &- t26 * (-t606 + t133) + s26 * (t54 + t634 - t727) + t275 + t302 -
     & t393 * s45) * s16 + s45 * (-t121 + t627) + (s56 * t786 - t153) * 
     &s46 - t149 + t252 + t794 + t303 - t632 * t366) * s12 - t578 * ((t1
     &02 * t763 + t765) * s16 + t258)
      t102 = t147 + t524
      t271 = ddreal(14) * s56
      t294 = ddreal(19) * s45
      t308 = t134 + t66
      t309 = t121 + t513
      t312 = s45 * t309
      t313 = ddreal(16) * s34
      t105 = -(s16 * (s45 - t105 - t565 - t791 - t806 + t144) - s26 * (t
     &271 + t805) + (-t675 + t780 - t907) * s46 - t116 * (t123 + t66) + 
     &t15 * t102 - t16 * t262 + t26 * (t560 - t148 + t176 + t231)) * s12
     & + s46 * ((-ddreal(21) * s26 + t632) * s45 - ddreal(11) * t86) + (
     &(-s45 * t815 - t313 - t451) * s26 - t738 + t133 * t2) * s56 - t15 
     &* (-t122 - t869 + t685 + t754 - t539) + t16 * (s26 * (t134 + t133 
     &+ t769) + (t153 + t553) * s46 - t149 + t248 - t692) + t26 * ((t148
     & + t905) * s34 + t124 + t252 + t299 + t613 + t687 + t688) + t61 * 
     &s46 * t298 - ddreal(13) * t312 + ddreal(9) * s26 * t889 - ddreal(9
     &) * (-t121 + t148 + t351) * s56 - (-(t16 * t495 + t345 + t359 - t5
     &97 + t697) * s16 - (-t667 - t158) * s34 - (-t257 - t154) * s45 + s
     &46 * (t526 + t803) - t26 * t308 - ddreal(14) * t121 - s25 * t26 * 
     &t273 - s26 * (-t294 + t393 + t271 - t679)) * s16 - t67
      t144 = t26 * t71
      t105 = -t19 * t510 * t671 + s12 * t105 + t144 * t261
      t154 = t751 * t166
      t166 = t15 - t725
      t211 = -s16 * t219 - s26 * t245 + t211 - t594 - t789
      t219 = t242 + t417
      t245 = t848 - t847
      t261 = t923 + s34
      t256 = s25 * (-s34 * t916 - t16 * t256 - t26 * (epinv * t637 + s56
     &))
      t273 = t751 * t211
      t298 = t147 * epinv
      t299 = ((-t237 - t848) * s56 - t153) * s34 + t247 - t688 + t946 + 
     &t208 * (-s45 + t397)
      t194 = (-t66 - t962) * s26 + t148 * t194 - t153 * t426 - t753 + (s
     &26 * (-t237 - t406) + (-t406 + t241 + t242) * s56 + t847 * t379) *
     & s25 + ((s45 + t397) * s56 + t718) * s26 - t969 - t217 * t148
      t194 = -t116 * ((s34 * t255 + t298) * s26 - s34 * (t134 + t298) + 
     &(-s46 * t949 + t134) * s56 + t121 * (s46 - t241)) - t15 * t405 * t
     &306 - t16 * t299 + t26 * t194 + t61 * ((t231 + t721) * s45 + (-t36
     &6 - t134 - t719 + t852) * s46 - t86 * t177) + ((-t219 * t61 - t26 
     &* (t1 * t246 + t405 + t848) + t161) * s16 - t116 * s45 * t261 + s3
     &4 * ((t2 - t959) * s56 - t632 - ddreal(7) * t847) + t16 * (t151 + 
     &t712 - t736) - t26 * ((t237 + t848) * s56 + t718) + t61 * (-t513 +
     & t415) + s26 * ((ddreal(9) - t817) * s34 - t15 * t245 - t803 + t54
     & * t27) + s46 * (t397 + t417) + t256) * s16 + s45 * ((-s56 + t237)
     & * s34 - t121 * t396 + t415) + (t712 + t858 + t721) * s46 + t148 *
     & (t406 - t241) + t749 + t273 + t231 * (s26 + t850 + t237) + ddreal
     &(7) * t147 * t320 + ddreal(9) * s34 * (t366 + t133) + ddreal(10) *
     & t86 * t940
      t217 = t578 * t510
      t255 = t513 + t133
      t299 = t116 * t1052
      t306 = t61 * t811
      t314 = ddreal(23) * s34
      t315 = t387 * s45
      t321 = t335 - t148
      t211 = s13 * t211
      t322 = (t392 - t153) * s56
      t326 = t133 * t128
      t328 = t86 * t74
      t333 = t121 - t148
      t336 = t392 - t351
      t338 = (-t153 + t66) * s46
      t339 = -t134 + t351
      t341 = ddreal(13) * s26
      t342 = t341 + t257
      t348 = -t578 * t510
      t353 = t16 * t694
      t349 = t349 * t121
      t355 = ddreal(12) * t366 * t74
      t210 = ((t116 * t871 + t26 * t416 - t61 * (t147 - t123 - t66 + t86
     &) + t184 + s16 * (t451 + t526 + t502 + t501) + t148 - t176 - t231 
     &- t560 - t353 + t291) * s12 + t116 * (-s56 * t339 + t869) + t16 * 
     &((-t148 + t627) * s26 + t122 - t686 - t776) + t26 * ((t121 + t210)
     & * s25 + (-t416 - t148) * s45 - s46 * t547 + t124 - t304 + t690) -
     & t61 * ((-t392 - t66 + t231) * s26 + t948) + ((t359 + t111 - t161 
     &- t306 + t501) * s16 + (-t212 - t588) * s46 + t61 * t66 - ddreal(1
     &1) * t102 + t276 + t604 + s26 * (t427 - t155 + t288 - t621) - t153
     & + t275 - t513) * s16 + t252 - t277 + t282 + t283 + t566 + t569 + 
     &ddreal(7) * t639 - ddreal(7) * t349 - t355 + t392 * (-t682 - t780)
     & + t598 * s34) * s12 - t578 * t784
      t357 = epinv * t162
      t358 = t116 * t875
      t360 = (-t406 - t848 + t405) * s45
      t364 = t747 * s56
      t369 = s45 * ((epinv + ddreal(1)) * s25 + s34)
      t372 = t2 * s26 * (t364 + t923)
      t374 = s45 + t62 + t190
      t375 = t2 * (-t151 + t121 - t351) - t61 * (-t66 + t148 + t133) - d
     &dreal(9) * t102 - t1111 + s26 * (-t15 * t811 - t727 + t804) + t203
      t378 = t524 * t199 - t122 - t643 - t665 + t686 + t693
      t380 = s56 * (t848 + t405)
      t382 = -t147 + t335
      t383 = (-t712 + t721) * s46 + t689 + t749 + t947 + t963
      t384 = (t1103 - t450) * s56 - t636 * t228
      t385 = ddreal(11) * s26 * (t147 * t922 + t86)
      t388 = ddreal(12) * t366 * (s45 * t922 - t405)
      t390 = (-s16 - s26 + s34 - s56 + t357 - t345) * s12 + t116 * ((-t4
     &06 + t177) * s26 - t740) + t15 * t650 + t16 * (s56 * (-t417 - t405
     &) - t123) - t26 * t369 + t61 * ((-s46 - t703) * s25 - t134 - t402 
     &- t721 + t887) + t1122 + s16 * ((t427 + t54) * epinv - t592 + ddre
     &al(7) * t657 - t805 + t358) - s25 * t445 + t154 - t720 - t360 + t3
     &72
      t391 = t217 * t801
      t141 = (t390 * s12 + s34 * t384 + t116 * ((t134 + t380 - t323) * s
     &45 + t776 - t151 * t387) - t15 * t405 * t382 - t16 * (epinv * t378
     & + (-t296 - t623) * s46 + t122 + t124 - t640 + t316 * t141) + t26 
     &* t383 - t61 * (s34 * ((-s25 - t228) * s34 - t134 + t335 - t237 * 
     &t246) + ((-t121 - t553) * epinv + t121) * s46 + t967 - t942 * t379
     &) + ddreal(7) * s34 * (-t524 + t316) + ddreal(7) * t121 * (t242 - 
     &t848 - t405) + ddreal(7) * t147 * (-s34 + t848) - ddreal(10) * t86
     & * (s34 + s45 - s46) + ((epinv * t374 - t212 + t597 + t621 - t631 
     &- t635) * s16 + (t313 + t682) * s56 - t16 * (t66 + t250) - t61 * t
     &773 - ddreal(12) * t121 + epinv * t375 + s26 * (s34 * t815 - t204 
     &- t451 + t681) + s56 * t563 - t153 - ddreal(10) * s34 * t924 + ddr
     &eal(10) * t147) * s16 + t153 * (t781 - t406 - t848) + t273 + t513 
     &* t744 + t207 * t320 + t415 * s45 + t385 + t388) * s12 + t391
      t250 = t848 + t527
      t296 = s34 - t418
      t383 = t406 * t2
      t384 = t15 * (s46 + t848)
      t390 = t207 * t1
      t391 = s46 + t237
      t393 = t736 - t415
      t394 = epinv * t811
      t398 = t116 * (s45 - t848)
      t403 = -t61 * t394
      t407 = (t111 + t501) * t1
      t409 = ddreal(23) * epinv
      t416 = t238 * t1
      t420 = t153 * t396
      t421 = ddreal(9) * s45
      t425 = t61 * t721
      t431 = t16 * t720
      t434 = t406 * t153
      t440 = (-t713 + t712) * s56
      t442 = t1 * t124
      t444 = (t405 + s26) * s56
      t298 = s34 * (t444 + t298)
      t445 = t242 * t121
      t446 = (t887 - t736) * s56
      t447 = -t415 * t811 + t445 + t446
      t448 = t177 * s45
      t449 = -epinv * t135 + (-t847 * t441 - t351) * s34 + s46 * (t717 -
     & t448) - t753 + t392 * (t417 + t177)
      t453 = s34 * t287 + t66 * (-t242 + t405)
      t454 = ddreal(12) * t715 * s45
      t298 = t116 * t453 - t16 * t449 + t26 * (((-s45 + t406) * s26 + (t
     &406 - t241 - t242) * s56 + t106 + t962) * s25 + (-t524 + t712 - t7
     &21) * s46 + t328 - t749 - t963 + t442) - t61 * ((t231 - t721) * s2
     &6 + s34 * (-t151 - t655) - t921 - t405 * t153 - t133 * t177) + (((
     &t730 - t16) * s34 + s46 - t398 - t403 + t407) * s16 + (-t521 - t42
     &1) * s46 + t26 * (s26 * (-t166 * s45 - t854 * s46 + t848 * t15) + 
     &t627) + ddreal(11) * t847 * t74 + ddreal(10) * s45 * (-t848 - t405
     &) + ddreal(10) * t736 + s25 * (-t26 * (s46 * t711 + s45 - t848) + 
     &t568 + t416) - t717 + t420 + t86 * (-ddreal(9) + t409) + t425 - t4
     &31 + t405 * (-t201 + t680) - t524 * t2) * s16 + (epinv * t321 + t1
     &48 + t153 + t316) * s34 + t134 * (s45 + t399) - t273 + t434 + t440
     & - ddreal(11) * t887 * s26 - ddreal(9) * t298 - ddreal(7) * t447 -
     & ddreal(14) * t402 * t74 + ddreal(10) * t405 * t336 - t454 + t86 *
     & (epinv * t342 - t427)
      t296 = ((-(-t26 * t912 - t250) * s12 + t116 * t393 - t16 * (s34 * 
     &t391 - t717) - t26 * (s26 * t296 - t209 - t390) - t61 * ((t242 + t
     &177) * s46 - t410 + t524 - t721) - t372 + ddreal(7) * t706 - s16 *
     & ((t451 - t588) * epinv + t26 * t296 + s45 - t383 - t384) + (-s34 
     &- t406 - t848 - t177) * s45 - t154 - t231 + t720 + ddreal(10) * t4
     &02) * s12 + t298) * s12 - t348 * t377
      t298 = -ddreal(10) * t716 + ddreal(10) * t736
      t377 = t228 * s46
      t447 = s34 * (t524 + t377)
      t449 = t415 * t199 + t242 * t66 + t135 + t776
      t245 = t116 * t449 - t15 * t133 * t728 - t16 * ((t209 - t720 - t71
     &8 - t415) * s26 + (-t316 + t1168) * s46 - t689 + t753 + t513 * t72
     &3 + t405 * t66) - t2 * t447 + t26 * ((t849 - t736) * s25 + (t390 +
     & t720 - t323 - t745) * s45 + s46 * t814 + t153 * (s25 - t417) - t4
     &42 + t231 * t245) + t61 * ((-t719 - t736) * s46 - t967 + t231 * (s
     &45 - t228)) - ddreal(13) * t328 + t385 + t388 + (((ddreal(10) - t4
     &08) * s34 - s46 + t398 + t403 - t407) * s16 - ddreal(11) * epinv *
     & t102 + (t405 * t26 + t204 + t427) * s34 - t425 + t431 - ddreal(13
     &) * t627 + s26 * ((ddreal(21) - t938) * s34 + (-t16 + t750) * s46 
     &- ddreal(12) * t848 + t54 * t166) + t256 - t420 + t717 - ddreal(9)
     & * s46 * t867 - t298) * s16 + (t335 - t720 + t962 + t448) * s34 + 
     &(t713 - t712) * s56 + t134 * (t399 + t237) + t273 - t434 + ddreal(
     &7) * (s46 * t957 - t852) * s34 + ddreal(7) * t445 + ddreal(7) * t4
     &46 + ddreal(9) * t133 * t387 + t313 * t366
      t138 = (((-t26 * t417 + s34 - t527 - t848) * s12 + (t138 + t526) *
     & s34 - t116 * (-t913 + t736) + t16 * (-t134 - t717 - t719) + t26 *
     & s25 * (t847 + t229 + t364) - t61 * (s45 * t638 + t209 + t402 + t4
     &10 + t721) + t372 + s16 * (-t26 * t920 + s45 - t383 - t384 - ddrea
     &l(7) * t918) - (s25 - t406 - t848 + t405) * s45 + t154 - t720 - t8
     &51) * s12 + t245) * s12 - t348 * t701
      t245 = (ddreal(1) + t732) * s34
      t256 = t133 * t563
      t323 = (-s26 * t343 - s56 * t379) * s25
      t128 = s34 * ((t341 + t257 + t803) * s26 - ddreal(10) * t351) + t1
     &16 * (s34 * (t153 - t852) + t66 * (s34 + t242) - t151 * t949) - t1
     &5 * ((t553 + t402) * s46 - t135) + t16 * (t148 * t903 + t640 - t75
     &3 + t904 + t948) + t26 * ((t286 - t285 + t284) * epinv + t149 + t3
     &38 - t442 + t692 + t323) - t61 * (s46 * t393 + t317 * t66 - t873 +
     & t231 * (-s34 + t228)) - ddreal(16) * t328 + t454 - (-((ddreal(12)
     & - t408) * s34 + t116 * (t242 - t848) - t407 - t417) * s16 - s34 *
     & ((ddreal(7) - t725) * s46 - t709) - t16 * (-t406 * t268 + t720) -
     & t160 * t26 + t298 + t425 - ddreal(9) * t887 - s25 * (-t26 * (t400
     & + t848 - t847) - t416 - t427 - t245) + t712 - t717 - t86 * (ddrea
     &l(25) - t959) + ddreal(16) * t352) * s16 - (t335 + t720 - t962 - t
     &448) * s34 + t273 - t434 - t440 + t415 * t128 + ddreal(7) * t446 +
     & ddreal(7) * t256 + ddreal(7) * t445 + ddreal(11) * s45 * (t392 * 
     &t922 + t134) + ddreal(9) * t524 * t735
      t128 = ((-(t26 * t638 - t250) * s12 - (t675 + t697) * s34 + t116 *
     & (t133 + t736) + t15 * t297 - t16 * s46 * (-t405 - s25 + t364) - t
     &26 * (s25 * (-s45 + t847) + (-t418 - t405) * s26 - t390) + t61 * (
     &t380 - t887 + t410) - t372 + s16 * (t15 * t944 + t406 * t2 + t26 *
     & t920 - t242 - t632) - t154 + t360 + t720 + t231 * t396) * s12 - t
     &128) * s12 + t348 * t787
      t297 = ddreal(9) * t392
      t298 = s26 * (t121 + t66)
      t127 = epinv * t221 + (-(t116 * t255 + t134 * t15 - t16 * (-s56 * 
     &t603 - t208 - t209) + t26 * (t581 + t951 + t66 + t207) - ddreal(7)
     & * t654 + s16 * (-ddreal(7) * t456 + t48 + t310) + t148 - t176 + t
     &231 + t297) * s12 + t116 * s34 * (t153 + t66) - t15 * ((t134 - t15
     &3) * s56 - t135 + t247 + t208 * t127) + t16 * (-s56 * t333 + t640)
     & + t26 * (-t298 + t149 + t323 + t692 + t754) - t61 * ((t121 - t66)
     & * s45 + t277 - t566 - t693) - ((t198 - t670 + t48) * s16 + t15 * 
     &s45 * t127 + (-t816 + t679 - t791) * s34 + t26 * t853 + t61 * t302
     & + t529 - t683 + s25 * (s34 + t427 + t238 + t195) + s26 * (t526 - 
     &t679 + t806 + t450)) * s16 - (t335 + t148) * s34 + t211 + ddreal(7
     &) * t121 * t456 + ddreal(7) * t256 - ddreal(7) * t869 + ddreal(11)
     & * (t392 + t134) * s45 + ddreal(11) * t366 * s34 + t86 * (-t679 + 
     &t806) + t513 * (t421 + t780)) * s12
      t71 = t26 * (t456 * t71 + t348 * (s16 * (s34 - s45 - s56 - s46) - 
     &s56 * t295 - t134 - t66 + t859 - t87 - t90)) + t127
      t87 = s25 * (t387 + t414)
      t90 = s34 * t373
      t127 = ((-s46 - t632 + t501 - t263) * s16 + t61 * (-s26 * (-t16 * 
     &t435 + t727) + t66) + t276 + ddreal(12) * t627 + t899 + t153 + t87
     & - ddreal(13) * t524) * s16
      t256 = t147 - t66
      t295 = t116 * t351 * t199
      t302 = t16 * (s34 * t308 - t122)
      t264 = (((t310 - t526 - t805) * s34 + t26 * t334 + t61 * (t324 + t
     &392 + t123) + t299 + s16 * (t249 + t526 - t588 + t501) - t176 + t3
     &16 + t337 - t458 + t648 + t366 * t2 + t340) * s12 + s34 * ((-t565 
     &+ t450) * s56 - ddreal(11) * t121) + t26 * ((-t147 - t153 + t553) 
     &* s46 + t124 + t605) - t61 * (s26 * t256 + t326) + (-t134 - t66 + 
     &t86 - t570) * s25 + s46 * (t121 - t883 + t90) + t252 + t748 - t969
     & + t127 + ddreal(12) * s26 * ((s34 - s56) * s45 + t134) - ddreal(7
     &) * s45 * t1052 + ddreal(7) * s56 * t1052 - t295 - t302) * s12 + t
     &264 * t578
      t308 = (t121 + t553) * s46
      t317 = t153 * t246
      t139 = ((-t116 * (-t151 + t581) - t61 * (s56 * t441 - t259 + t649)
     & - t184 + t353 - s16 * (t427 + t54 + t168) - t148 + t176 + t560 + 
     &t111 * s45 - t260) * s12 + t15 * s34 * t382 + t116 * s45 * ((-s34 
     &+ s56) * s56 + t335) + t16 * t378 - t26 * (s46 * (t153 - t66 - t85
     &) + t566) - t61 * (s34 * (s45 * t246 + t86) + t248 + t613 + t308) 
     &+ t355 + (-s16 * t374 - t375) * s16 - t252 + t322 - t689 - t868 - 
     &t317 + ddreal(7) * t349 + ddreal(7) * t612 + t392 * (t682 + t588))
     & * s12 + t139 * t510
      t259 = s34 + t357
      t323 = -s26 * t935
      t298 = (-t121 + t351 + t524 - t436 + t415) * s25 + (t153 + t351) *
     & s34 - t298 + t749 + t963 + t721 * s46
      t334 = t851 * t637
      t337 = ((s56 - t417 - t177) * s46 - t397 * t441) * s26 + s45 * (-t
     &66 + t812) - t135 - t566 + t775 - t334
      t305 = epinv * ((s45 * t305 - t134) * s46 + (t208 - t134) * s56)
      t341 = (-t153 - t852) * s25
      t343 = (t121 + t134 - t720) * s56 + t809 + t341 - t852 * t857
      t195 = s26 * (((ddreal(12) * epinv - ddreal(14)) * s34 - t682) * s
     &46 + ddreal(17) * t719) + t116 * ((-t958 - t66 - t133) * s34 + t20
     &8 * t177) + t15 * (((-s46 + t228) * s25 + t230 + t862) * s34 + t53
     &8 - t748 + t147 * t426) + t16 * t343 + t2 * t305 - t26 * t298 + t6
     &1 * t337 + (t153 * t251 - t868) * epinv + ((-t26 * t711 * t246 + t
     &116 * t219 + t16 * t940 - t26 * t886) * s16 - t116 * t237 * t533 +
     & s34 * ((ddreal(14) * epinv - ddreal(13)) * s46 + t727 + t397 * t6
     &96) - t15 * t712 + t26 * (t853 - t415) - t61 * ((-s56 + t417) * s4
     &6 + t230) - ddreal(9) * t524 - ddreal(9) * t887 + s25 * (t15 * t41
     &9 - t223 * t61 + s34 + t195 + t323 + t937) + ((-epinv + t2) * s45 
     &+ (t15 - t938) * s56 + ddreal(12) * t419 + ddreal(9) * t918 + t323
     &) * s26) * s16 - t273 + t685 + ddreal(7) * t121 * t728 + ddreal(7)
     & * t869 - ddreal(7) * t133 * (s25 - s46 + t397) + ddreal(7) * t406
     & * t134 - ddreal(10) * t86 * (s56 - t228) - ddreal(9) * t513 * t74
      t48 = (((t259 - t345) * s12 + t116 * (-t228 * t227 - t513) + t16 *
     & (s26 * t619 + t281 - t351) + t2 * t228 * t199 - t26 * ((s34 + t39
     &9 - t417) * s45 + t207 + t581 + t66) + t61 * (epinv * t173 - t564)
     & - t297 - s16 * (-t16 * t741 - t358 + t405 + ddreal(7) * t417 + t4
     &8) - s34 * t731 + t154 - t720 + ddreal(7) * t913 - t719 * t15) * s
     &12 - t195) * s12 + t217 * t593
      t195 = t335 - t524
      t219 = s45 - t223
      t227 = ddreal(32)
      t246 = epinv * t227
      t281 = -ddreal(27)
      t297 = t392 * t74
      t245 = ((t26 * t930 * s25 + (t696 - t408) * s34 + (t61 + t710) * s
     &45 - t636 + ddreal(10) * t968 + ddreal(13) * t397 - ddreal(14) * t
     &406) * s16 + s34 * (-ddreal(17) * s34 + ddreal(23) * s46) + ((ddre
     &al(21) - t970) * s46 + t142 - t201) * s45 + ((-ddreal(19) * epinv 
     &+ ddreal(11)) * s45 + s46 * t289 + t405 * t281 + t730 * s56) * s56
     & - t15 * (t66 - t720 + t962 + t718) + t26 * (t121 * (-t2 + t938) +
     & t415) - ddreal(24) * s34 * t968 + s25 * (-t15 * t219 + t16 * (s26
     & * t930 - t245) - s56 - ddreal(10) * t391) + s26 * ((-ddreal(21) +
     & t246) * s56 - ddreal(19) * t417 + (ddreal(31) - t409) * s34 + (dd
     &real(9) - ddreal(13) * epinv) * s45)) * s16 + s56 * (-s25 * t662 +
     & t230)
      t259 = (-t259 * t26 + t212 + t363) * s12 + t116 * t66 - t15 * (-t1
     &77 * t657 - t134 + t887) + t16 * t369 + t26 * ((t726 + t846 + t405
     &) * s45 - t154 + t231 + t720) + t389 * t61 + ddreal(17) * t713 + d
     &dreal(12) * t746 - ddreal(9) * s56 * t261 - ddreal(10) * epinv * t
     &195 + ddreal(11) * s26 * t735 + ddreal(7) * t209 + ddreal(7) * t12
     &1 - ddreal(7) * t230 - s16 * (s16 * t931 - (-t730 + t61) * s25 - (
     &t26 - t725) * s45 + ddreal(12) * t292 - t405 - ddreal(14) * t417 -
     & t658 - ddreal(11) * t673) + t228 * (-ddreal(18) * s56 + t780) - t
     &572 * s34
      t261 = t348 * (t952 + t121 - t917 + t437 + t876 + t954 - t953)
      t177 = (t259 * s12 - (((ddreal(25) - t246) * s26 + t572) * s34 - t
     &366 * t116 - ddreal(26) * t919) * s56 - t15 * (-t812 * s45 - t335 
     &* t219 - t122 + t135 + t334 + t566 + t775) + t16 * ((t335 - t153) 
     &* s34 - t749 - t947 - t963 + t712 * s46) - t2 * (-t341 - t795 + t5
     &36) - t26 * ((s45 * t731 + t231) * s56 + t124 * t27 + t273 - t665 
     &+ t237 * (t641 + t134) - t847 * t66 + t228 * t153) + t61 * (t513 *
     & (-s56 + t406) - t748) + t696 * t237 * t309 - ddreal(14) * s45 * (
     &t231 + t444) + ddreal(12) * s34 * ((t177 + s45) * s26 + (s46 * t92
     &2 - t397) * s34) + ddreal(10) * t228 * (-t121 - t134 - t66 - t661)
     & + ddreal(10) * t308 - ddreal(10) * t868 - ddreal(24) * t297 - ddr
     &eal(16) * t147 * (s56 - t241) - ddreal(16) * t543 - t245 + ddreal(
     &23) * t377 * t74 + ddreal(11) * s56 * ((t929 + s34) * s34 + t121 +
     & t861) + ddreal(19) * t852 * t387 - ddreal(9) * t866 - ddreal(9) *
     & t776 + ddreal(7) * t148 * (t364 + t229) + t677 * t134) * s12 + t2
     &61
      t219 = -t153 - t524
      t229 = t300 + t133
      t89 = (-t26 * t419 - t250) * s12 + t116 * epinv * t1052 - t15 * (t
     &392 + t719) + t16 * (-t366 + t718) + t26 * ((-s26 - s56 - t417) * 
     &s25 + t694 + t147 * t396) + t61 * (t151 - t121 + t230 + t862 + t37
     &7) - (-(t526 - t588 + t501) * epinv - s34 + t394 + t412 - t413 + t
     &505 - t856) * s16 - (-s34 - t89) * s25 + t154 - t280 + t344 - t66 
     &- t720 - t740 + t228 * (-t588 + t635)
      t84 = (t89 * s12 + t116 * (t366 * s16 + t613 + t693 + t724) + t15 
     &* s46 * (t208 + t900) + t16 * ((-t335 + t66 + t148) * s16 + t936) 
     &+ t26 * (s16 * t229 + s26 * t229 + (t209 + t66 + t148 + t207) * s4
     &5 - t121 * t646 - t123 * t810 + t149 + t304 + t67) + t61 * ((s16 *
     & t84 + t366) * s26 + t770 + t231 * t107) + epinv * (-t26 * ((t147 
     &+ t231 + t905) * s46 - t124 + t640) - t326 * t61 - t295 - t302 - (
     &t324 + t134) * s25 - (-(s46 + ddreal(7) * t435 - t621) * s26 + t61
     & * t256 - t153 + ddreal(12) * t339 - t87 - t293 * (t597 - t675)) *
     & s26 - s46 * (t883 - t90) + t127 + t688 - t969 + ddreal(7) * t134 
     &* t435 + t795 * t2) - s16 * t219 - t219 * s26 + s45 * t123 - (t315
     & + t134) * s45 - ((-s16 - s26 + s34 + s45 + s56) * s34 - t123) * s
     &46 + t273 - t303 - t539 - t569 + t513 * (t697 + t114)) * s12 + t33
     &1 * t578
      t67 = (s26 * (t636 + t103) + (-s34 * t227 - ddreal(26) * s45 + t62
     &) * s56) * s26 + t15 * s25 * (-s34 * t233 + (s34 - s56 - s46) * s4
     &5 + t234) + t16 * (s46 * (t153 - t85) + t566) + t26 * (t124 + t252
     & + t868 + t689 + t688 + t317 - t754) - t312 * t696 - ddreal(23) * 
     &t297 - ddreal(12) * (t335 - t650) * s34 + ddreal(10) * s26 * ((-s2
     &5 + s34) * s45 + t121 + t134 + t66) - ddreal(11) * s56 * ((s34 + s
     &45) * s56 + t151) - ddreal(7) * t148 * t320 + ((t54 - t161 + t565 
     &+ t667 + t158 + t103) * s16 + s26 * (s56 * t227 + ddreal(19) * s46
     & - t314 - t681) + (s34 * t281 + t212 - t294 + t791) * s56 - t15 * 
     &(t148 + t789) + t26 * (t134 + t203) - t815 * t287 - ddreal(24) * t
     &151) * s16 + t67 - ddreal(19) * t121 * t387 + t133 * (t808 + t271)
      t67 = t144 * t162 + ((s26 * (ddreal(17) * s46 + t157) + (-t806 - t
     &501) * s34 + (t874 + t636) * s56 + t15 * (-t147 - t172) - t26 * (t
     &560 - t148 + t176) - t587 + ddreal(7) * t66 + s16 * (s16 - s34 + t
     &218 + t1064 + t798 + t772) + ddreal(10) * t195) * s12 + t67) * s12
     & - t439 * t578
      t87 = s12 - s16 - s26 - s34 - s45 - t200
      t89 = (s12 + t367) * s12 + t348
      t90 = s12 + s34 + s45 + s46
      t7 = t20 * t50 + t7 * ((s12 + t726 + t137 + t223 - t345 - t589) * 
     &s12 + t578 * t985)
      t7 = t165 * (t18 * (-t554 * t48 - t559 * t82 - t571 * t210 + t574 
     &* t139 - t633 * t141 - t700 * (((s12 * t638 + t16 * ((s26 - t847) 
     &* s34 + (-s46 - t406) * s56 - t134) + t26 * ((epinv * t301 - s25) 
     &* s26 + (t237 - t418) * s56 - t721 - t736 + t745) + t61 * ((-s25 +
     & s34 + t397) * s46 - t560 + t802) - s16 * (-t16 * t638 + t242 * t6
     &1 + t26 * (t406 + t418 + t848) + t405) - t316 - t720 - t865 + t154
     & - t392 * t166) * s12 + t194) * s12 + t217 * t882) - t926 * t138) 
     &+ t8 * (-t139 * t551 + t210 * (-t910 + t768))) + t18 * t7 + t319 *
     & (t8 * (t428 * ((s12 + t365 - t54) * s12 - t332 * t578) + (t356 - 
     &t386) * (s12 * (-t41 + t345 + t142) + t19 - t361)) + t34 * t21 * t
     &90) + s34 * (t89 * t983 + t902 * (t89 * epinv + t319 * t90)) * t8
      t2 = t622 * ((((t212 - t670) * s56 - t116 * t232 + t173 * t61 - t2
     &6 * t966 + s16 * (-s34 + t168 + t310 + t451) - t133 + t148 - t176 
     &+ t184 - ddreal(7) * t151 + t260) * s12 + t116 * t208 * t129 - t15
     & * ((t516 + t524) * s56 + t685) + t16 * (t121 * t249 + t122 - t536
     & - t686) - t2 * ((-t793 - t134) * s46 + (t208 - t134) * s56) + t26
     & * ((-t66 - t85) * s46 + t566) + t61 * (t248 + t247 - t693 + t739 
     &+ t613) - (-(t198 - t293 + t190) * s16 - (-t218 - t204) * s34 - t1
     &16 * t534 + t15 * t153 - t61 * (t66 - t148) - t1083 - t604 + t683 
     &- s26 * (s45 + t636 - t201 + t804) - t203) * s16 - s45 * (s45 * t2
     &51 - t134) + t252 + ddreal(7) * s56 * t774 + ddreal(7) * t642 + t8
     &6 * (-t257 - t697 - t806)) * s12 - t578 * (s16 * (-t785 * s16 - t2
     &6 * (-(t199)**2 + t153) + s34 * (s34 - s45 - t764) + t238 * t763) 
     &+ t258)) * t165
      t2 = t18 * (t1012 * (-t136 * t484 - t278 * (epinv * t136 + s12 * (
     &s25 * (t320 - t345 - t589) + t113 - t325 + t118 + t117) - t578 * (
     &s25 * (s26 - t200 - t589) + t109 - t436 - t438 - t853 - t112) - t1
     &32) + t307 * (-epinv * t65 + (-s25 * (-s46 + t423 - t589) - t311 *
     & s46 - t101 + t119) * s12 + t578 * (s25 * (-t199 - t589 + t60) + t
     &101 - t106)) + t318 * t91 + t452 * (epinv * t91 - t26 * (t170 - t1
     &31)) - t485 * t65 - t2 - t28 * t50) + t73 * ((s16 + s26 + s56 - s4
     &6 - t75 + t39) * s12 - t578 * t76)) + t8 * (-t38 * t465 + t69 * ((
     &s12 + t397 + t726 + t137 - t400 - t850) * s12 - t120 * t578) + t26
     &6 * t995 * t139 * t1012) + t56 * (t362 * (s12 * t78 + t25 * (t70 -
     & t19)) + t422 * (((-t934 + t225 - t200) * s12 + epinv * t78 + t169
     & * t83) * s12 - t350 * t578))
      t15 = t578 * t1012
      t20 = t508 * t264 - t558 * (((-(t363 - t1099) * s12 - t16 * t255 -
     & t26 * t389 - t61 * (t649 + t123 + t66) - t184 + ddreal(7) * t650 
     &- s16 * (t111 + t451 + t526 - t588 + t501) - (s25 - s26 - s56) * s
     &45 - t148 + t176 - t299 + ddreal(10) * t86) * s12 + s26 * (s34 * t
     &342 - t290) + t116 * t66 * t74 - t16 * (t333 * s26 - t513 * t373 -
     & t135 + t247) + t26 * ((-t366 - t121 + t153 + t351 + t524) * s25 -
     & t124 + t149 + t692 + t338) + t61 * (-s26 * t863 + t326) + (-(t359
     & + t111 - t306 - t632 + t501) * s16 + (t313 - t412) * s56 - ddreal
     &(10) * t121 - ddreal(12) * t134 - ddreal(10) * t315 + t504 + s25 *
     & (-t143 * t16 + s34 - t43) + s26 * (-t427 + t155 - t288 + t314) + 
     &t153 + t513 + ddreal(11) * t516) * s16 + s34 * t321 + t566 - t687 
     &+ t211 + t322 - ddreal(14) * t328 + ddreal(10) * s34 * t336 + ddre
     &al(7) * s46 * t339 - ddreal(7) * t759 + ddreal(7) * t868 - ddreal(
     &9) * s34 * t102) * s12 - t348 * t651) - t660 * t296 + t845 * t84
      t20 = t1012 * (t10 * t20 + t18 * (t443 * t71 - t561 * t128 + t575 
     &* t221 + t616 * (((-t16 * t567 + t26 * t262 + t61 * (-t147 - t655 
     &+ t123) + t184 + t684 - s16 * (s45 - t634 - t502 - t501) + t148 - 
     &t176 - t231 - t560 + t291) * s12 + s26 * ((-t427 - t658) * s34 + t
     &290) + t116 * (s45 * t901 - t543) + t16 * ((-t148 + t133) * s56 + 
     &t122 - t693) + t26 * (t124 + t286 - t285 + t284) - t61 * (-s46 * t
     &1052 + t531 * t66 - t539 + t685) + ((s46 + t111 - t161 + t501 - t2
     &63) * s16 + (t412 - t780) * s56 - t16 * (-s26 * t268 + t148) + t26
     & * t627 - ddreal(9) * t147 - t267 - ddreal(9) * t86 + t153 - t513 
     &+ t275 + t276) * s16 + t252 + t566 + t569 - t277 + t282 + t283 - d
     &dreal(7) * s56 * t287 - ddreal(7) * t775) * s12 - t145 * t510) + t
     &733 * t177 - t799 * t67) + t669 * t8 * t105) * t165
      t9 = t279 * (t1011 * (t1012 * (t18 * (-t11 * t14 + t429 * (epinv *
     & t152 + s12 * (s25 * (t320 - t62 - t432) - t325 + t113 + t118 + t1
     &17) - t578 * (-t438 + t433 + t109 - t853 - t436 - t112) - t132) - 
     &t186 * t50) + t9 * t8 * t11) + t578 * (t1012 * (t10 * t240 * (epin
     &v * t87 - s12 + s34 + s45 + t164 + t243) - t18 * t236 * (-t26 * (s
     &26 - s34 - s45 - t223) + epinv * t202 + t235 - t239)) + t50 * (t37
     &6 + t115)) + t20) + t1012 * (t10 * (-t68 + t174 + t475) + t21 * (-
     &s34 * t464 + t50 * t81 + t53)) * t319)
      t14 = t1026 * (-t178 * t38 + t29 * t30 + t347 * ((t556 * t61 + t10
     &71 + (t520 - t589 + t60 + s16) * s16 + (t520 - t589 + s26) * s26 -
     & t147 - t63 + t340) * s12 - t354 * t578) + t395 * (s12 * (s12 * t1
     &07 + (-s46 + t49) * s46 + (t49 - t345 - s56) * s56 - t844) + t52 *
     & t346) + t404 * (((-t26 * t329 + t424) * s12 + t148 * t26 + t16 * 
     &t160 - t61 * (-s46 * (s45 + t405) + (-s46 - t242) * s56) - (t156 +
     & t591 - t159 + t220) * s16 + s26 * (-t156 + t159) + s46 * t941 - (
     &t850 + t847) * s56 + t66 - t736) * s12 - t411 * t578))
      t5 = t401 * t1031 * (t8 * (t1015 * (t192 * (s12 * (t26 * t363 - t9
     &24) - t19) + t216 * (s12 * (t107 * t927 - t26 * t728 + t241 + t242
     & - t412 + t847) + t5) + t22 * (t181 - t40) + t32 * ((s12 * t927 - 
     &t224 * t61 + t432 + t107 * (-t16 + t732) - t848) * s12 + t5)) - t1
     &25 * t699) + t14 + (-t197 * t42 - t222 * (epinv * t42 - t57) + t24
     &4 * ((t52 * t319 + t110 + t119) * s12 - t131) - t33 * (epinv * t11
     & - t57) + t370 * t152) * t1015 * t18 + t1015 * (t50 * (t59 * s34 -
     & t1016 * t12) * t18 + t699 * (t1 * (t1003 * t150 + t50 * (-t1001 +
     & t1002)) - t50 * t1020 * t27) * t21) * t1013)
      t4 = t5 + t9 + t1031 * t1030 * (t1012 * ((t175 * (t39 * t699 + t19
     &3) - t35 * ((s12 + t397 + t726 + t137 - t417 - t75) * s12 - t19)) 
     &* s34 + t10 * (t29 * t31 - t38 * t992) * epinv - t1024 * t8 * t105
     & * t165 + (-t150 * t807 - t972 * ((t100 - t52) * s25 + t171 - t626
     & + t758)) * t1013 * t699 * t21) + t50 * (t1012 * t1013 * (t12 * t9
     &75 + t253) * t18 + t1013 * (t21 * (t1012 * (t979 - t980 - t981 - t
     &982) + t188 + t191 - t64) + epinv2 * t4 * t1012) * t699 + (t1012 *
     & ((t993 - t1007) * epinv + t180 - t986) + t146) * t578 * t8 + t140
     & * t56) - t15 * (-t18 * t215 + t214 * t8) * t87)
      t5 = t1022 * propW34 * (t1023 * (t8 * (t274 * t55 * (s12 - s45 - s
     &46 - s56) * t319 + t270 - t72) - t96 + t15 * t1030 * (t21 * (-t100
     &5 * t319 + t430 * (-epinv * t319 + s12 + s16 + s26 - s34 - s45)) +
     & t319 * (-t226 + t1004) * t8) * s34 + t254 * (t1012 * (t269 * (-ep
     &inv * t98 - t131 + t170) - t483 * t98) + t466 * (-(t36 + t39) * s1
     &2 + t19))) + t92 * t1012 * (s16 * (t1 * (t998 - t1009) + t206 + t9
     &73 - t974 + t994) - t183 - t213) * t50 * t699 * t185)
      t2 = -t126 * t108 * t699 * t50 * t1022 * t1012 * t88 + t196 * (t10
     &31 * (t18 * t95 + t401 * (t265 - t51) * t699 * t21 + t272 * t167 *
     & t319 + t254 * t699 * t50 * t1012 * t839 * t1013 * s34 + t1030 * t
     &2 - t327 * (t486 * t99 + t58 * t97)) + t279 * (t1012 * (t1011 * t7
     & + t130 * t165 * (t58 * t991 + t99 * t990)) + t37 * t50 * t1011 * 
     &t1013 * t699)) + t104 * t4 + t5 / ddreal(3)
      t4 = t1031 * t196 * (-t26 * (t21 * (-t1 * (t633 + t733 + t926 + t3
     &07 + t278) + t484 + t485 + t571 + t574 + t799) + t8 * (t1 * (t995 
     &+ t996 + t997 + t964 + t800) - t330 - t386 - t551 - t669 - t768)) 
     &- t10 * (t1 * (t660 + t404) - t347 - t558) - t18 * (-t370 - t616) 
     &- t1 * (t561 + t429) * t21) * struc11Step4 * t1012
      result = mytree * (t196 * (t205 + t17 + t1031 * (t18 * t94 + t1030
     & * (t10 * (-t79 * (s12 * t13 - t80) + t46 * t467) - t376 * (-t77 +
     & t46)) + t381 * t1014 * t1032 * t24 + t93) + t1031 * (t1012 * (t10
     &11 * (t1015 * (t8 * (epinv * (t179 * t992 + t182 * t993) + t46 * (
     &epinv * t1007 + t986)) - t1024 * t10 * t611 * t165) + t1025 * t462
     &) + t1026 * (t1 * t474 * t988 + t1021 * t24) * t1014) + t1030 * (t
     &189 * (t1 * (t1027 + t977) + t1012 * (-t979 + t980 + t981 + t982) 
     &+ t1028) + t1029 * t47 * s34 + t1012 * t1017 * (t368 * t44 + t46 *
     & t976) - t1 * t371 * t1003 * t1012) * t1013 * t18)) + t126 * t23 -
     & t163 * t45 + ddreal(4) / ddreal(3) * t108 * t187 * t1022 * t1012 
     &* t88) + t3 * struc4Step5 + t6 * struc8Step5 + t2 * struc9Step5 + 
     &t4

           ampNonresonantLightFullMM = result
       end function ampNonresonantLightFullMM

       function ampNonresonantLightFullMP()
           implicit none
           type(dd_complex) :: ampNonresonantLightFullMP

           type(dd_complex) :: result

      result = ddreal(0)

           ampNonresonantLightFullMP = result
       end function ampNonresonantLightFullMP

       function ampNonresonantLightFullPM()
           implicit none
           type(dd_complex) :: ampNonresonantLightFullPM

           type(dd_complex) :: result

      result = ddreal(0)

           ampNonresonantLightFullPM = result
       end function ampNonresonantLightFullPM

       function ampNonresonantLightFullPP()
           implicit none
           type(dd_complex) :: ampNonresonantLightFullPP
           type(dd_complex) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           type(dd_complex) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           type(dd_complex) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           type(dd_complex) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           type(dd_complex) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           type(dd_complex) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           type(dd_complex) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           type(dd_complex) ::  t175,t176,t177,t178,t179,t18,t180,t181,t182,t183,t184,t185
           type(dd_complex) ::  t186,t187,t188,t189,t19,t190,t191,t192,t193,t194,t195,t196
           type(dd_complex) ::  t197,t198,t199,t2,t20,t200,t201,t202,t203,t204,t205,t206
           type(dd_complex) ::  t207,t208,t209,t21,t210,t211,t212,t213,t214,t215,t216,t217
           type(dd_complex) ::  t218,t219,t22,t220,t221,t222,t223,t224,t225,t226,t227,t228
           type(dd_complex) ::  t229,t23,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239
           type(dd_complex) ::  t24,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249,t25
           type(dd_complex) ::  t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t26,t260
           type(dd_complex) ::  t261,t262,t263,t264,t265,t266,t267,t268,t269,t27,t270,t271
           type(dd_complex) ::  t272,t273,t274,t275,t276,t277,t278,t279,t28,t280,t281,t282
           type(dd_complex) ::  t283,t284,t285,t286,t287,t288,t289,t29,t290,t291,t292,t293
           type(dd_complex) ::  t294,t295,t296,t297,t298,t299,t3,t30,t300,t301,t302,t303
           type(dd_complex) ::  t304,t305,t306,t307,t308,t309,t31,t310,t311,t312,t313,t314
           type(dd_complex) ::  t315,t316,t317,t318,t319,t32,t320,t321,t322,t323,t324,t325
           type(dd_complex) ::  t326,t327,t328,t329,t33,t330,t331,t332,t333,t334,t335,t336
           type(dd_complex) ::  t337,t338,t339,t34,t340,t341,t342,t343,t344,t345,t346,t347
           type(dd_complex) ::  t348,t349,t35,t350,t351,t352,t353,t354,t355,t356,t357,t358
           type(dd_complex) ::  t359,t36,t360,t361,t362,t363,t364,t365,t366,t367,t368,t369
           type(dd_complex) ::  t37,t370,t371,t372,t373,t374,t375,t376,t377,t378,t379,t38
           type(dd_complex) ::  t380,t381,t382,t383,t384,t385,t386,t387,t388,t389,t39,t390
           type(dd_complex) ::  t391,t392,t393,t394,t395,t396,t397,t398,t399,t4,t40,t400
           type(dd_complex) ::  t401,t402,t403,t404,t405,t406,t407,t408,t409,t41,t410,t411
           type(dd_complex) ::  t412,t413,t414,t415,t416,t417,t418,t419,t42,t420,t421,t422
           type(dd_complex) ::  t423,t424,t425,t426,t427,t428,t429,t43,t430,t431,t432,t433
           type(dd_complex) ::  t434,t435,t436,t437,t438,t439,t44,t440,t441,t442,t443,t444
           type(dd_complex) ::  t445,t446,t447,t448,t449,t45,t450,t451,t452,t453,t454,t455
           type(dd_complex) ::  t456,t457,t458,t459,t46,t460,t461,t462,t463,t464,t465,t466
           type(dd_complex) ::  t467,t468,t469,t47,t470,t471,t472,t473,t474,t475,t476,t477
           type(dd_complex) ::  t478,t479,t48,t480,t481,t482,t483,t484,t485,t486,t487,t488
           type(dd_complex) ::  t489,t49,t490,t491,t492,t493,t494,t495,t496,t497,t498,t499
           type(dd_complex) ::  t5,t50,t500,t501,t502,t503,t504,t505,t506,t507,t508,t509
           type(dd_complex) ::  t51,t510,t511,t512,t513,t514,t515,t516,t517,t518,t519,t52
           type(dd_complex) ::  t520,t521,t522,t523,t524,t525,t526,t527,t528,t529,t53,t530
           type(dd_complex) ::  t531,t532,t533,t534,t535,t536,t537,t538,t539,t54,t540,t541
           type(dd_complex) ::  t542,t543,t544,t545,t546,t547,t548,t549,t55,t550,t551,t552
           type(dd_complex) ::  t553,t554,t555,t556,t557,t558,t559,t56,t560,t561,t562,t563
           type(dd_complex) ::  t564,t565,t566,t567,t568,t569,t57,t570,t571,t572,t573,t574
           type(dd_complex) ::  t575,t576,t577,t578,t579,t58,t580,t581,t582,t583,t584,t585
           type(dd_complex) ::  t586,t587,t588,t589,t59,t590,t591,t592,t593,t594,t595,t596
           type(dd_complex) ::  t597,t598,t599,t6,t60,t600,t601,t602,t603,t604,t605,t606
           type(dd_complex) ::  t607,t608,t609,t61,t610,t611,t612,t613,t614,t615,t616,t617
           type(dd_complex) ::  t618,t619,t62,t620,t621,t622,t623,t624,t625,t626,t627,t628
           type(dd_complex) ::  t629,t63,t630,t631,t632,t633,t634,t635,t636,t637,t638,t639
           type(dd_complex) ::  t64,t640,t641,t642,t643,t644,t645,t646,t647,t648,t649,t65
           type(dd_complex) ::  t650,t651,t652,t653,t654,t655,t656,t657,t658,t659,t66,t660
           type(dd_complex) ::  t661,t662,t663,t664,t665,t666,t667,t67,t68,t69,t7,t70
           type(dd_complex) ::  t71,t72,t73,t74,t75,t76,t77,t78,t79,t8,t80,t81
           type(dd_complex) ::  t82,t83,t84,t85,t86,t87,t88,t89,t9,t90,t91,t92
           type(dd_complex) ::  t93,t94,t95,t96,t97,t98,t99

           type(dd_complex) :: result

      t1 = intHL0s25s26s34s56x1120D2eps0()
      t2 = ddreal(3)
      t3 = ddreal(2)
      t4 = t2 * s25
      t5 = t4 * propZ25
      t6 = -t3 + t5
      t7 = (gb)**2
      t8 = (gw)**2
      t9 = propZ25 * s25
      t10 = t9 * t7
      t11 = t8 * t6
      t12 = t11 + t10
      t13 = intHL0s25s26s34s56x1130D4eps0()
      t14 = intHL0s25s26s34s56x1210D2eps0()
      t15 = intHL0s25s26s34s56x1220D4eps0()
      t16 = intHL0s25s26s34s56x1310D4eps0()
      t17 = intHL0s25s26s34s56x1130D4eps1()
      t18 = ddreal(1) - epinv
      t19 = intHL0s25s26s34s56x1220D4eps1()
      t20 = intHL0s25s26s34s56x1310D4eps1()
      t21 = intHL0s25s26s34s56x1120D2eps1()
      t22 = intHL0s25s26s34s56x1210D2eps1()
      t23 = intHs16s25s26s34s56x1112D2eps0()
      t24 = ddreal(4)
      t25 = t24 - t5
      t26 = t8 * t25 + t10
      t27 = s25 + s26 + s56
      t28 = intHs16s25s26s34s56x1113D4eps1()
      t29 = intHs16s25s26s34s56x1114D6eps0()
      t30 = intHs16s25s26s34s56x1123D6eps0()
      t31 = intHs16s25s26s34s56x1213D6eps0()
      t32 = intHs16s25s26s34s56x1114D6eps1()
      t33 = intHs16s25s26s34s56x1123D6eps1()
      t34 = intHs16s25s26s34s56x1213D6eps1()
      t35 = intHs16s25s26s34s56x1112D2eps1()
      t36 = t32 * t18
      t37 = t36 - t29
      t38 = t26 * t37
      t39 = t17 + t19 + t20
      t40 = t18 * t39 - t13 - t15 - t16
      t41 = (t18 * (t33 + t34) - t28 - t30 - t31) * t26
      t42 = (t21 + t22) * t18
      t43 = t1 + t14 - t42
      t44 = epinv * t35
      t45 = (t44 + t23) * t26
      t46 = ddreal(6)
      t47 = intHLs160000x0111D0eps0()
      t48 = intHLs160000x0112D2eps0()
      t49 = intHLs160000x0121D2eps0()
      t50 = intHLs16s25s26s34s56x1112D4eps0()
      t51 = intHLs16s25s26s34s56x1121D4eps0()
      t52 = intHLs16s25s26s34s56x1141D6eps1()
      t53 = s46 + s56 + s45 - s12
      t54 = intHLs16s25s26s34s56x1231D6eps1()
      t55 = intHLs16s25s26s34s56x1321D6eps1()
      t56 = intHLs16s25s26s34s56x1112D2eps1()
      t57 = t2 * epinv
      t58 = s25 + s26
      t59 = intHL0s25s260s56x1012D2eps0()
      t60 = -s12 + s25 + s45 + s46
      t61 = intHL0s25s260s56x1013D4eps0()
      t62 = intHL0s25s260s56x1021D2eps0()
      t63 = intHL0s25s260s56x1021D2eps1()
      t64 = epinv * t60
      t65 = intHs160s26s34s56x1020D2eps0()
      t66 = s12 + s13 + s16 + s26 - s34 - s45 - s46
      t67 = s16 + s26 - s34 + s56
      t68 = intHs16s25s26s34s56x1110D2eps1()
      t69 = -t3 + epinv
      t70 = intHLs16s25s26s34s56x1131D4eps1()
      t71 = -t3 + t57
      t72 = s46 + s45
      t73 = t2 * t72
      t74 = t24 * s56
      t75 = s45 + s46 + s56
      t76 = t71 * s12
      t77 = t3 * t75
      t78 = intHs16s25s26s34s56x1110D2eps0()
      t79 = -s12 + s13
      t80 = s12 + s16
      t81 = t3 * s13
      t82 = t80 - t81
      t83 = s12 - t81
      t84 = t3 * s25
      t85 = s16 - t84
      t86 = t46 * s25
      t87 = (s16)**2
      t88 = s16 * t87
      t89 = s13 * t85
      t90 = s25 * t58
      t91 = s34 + s26
      t92 = t84 * t58
      t93 = s12 * s25
      t94 = s16 * s46
      t95 = ddreal(7)
      t96 = ddreal(5)
      t97 = (s25)**2
      t98 = s25 * t97
      t99 = t24 * s16
      t100 = t95 * s25
      t101 = t3 * s16
      t102 = t4 - t101
      t103 = s16 - t79
      t104 = (s34)**2
      t105 = t83 * t97
      t106 = s16 * t104
      t107 = t87 * t103
      t108 = s16 * s25
      t109 = t108 * t82
      t110 = (s26)**2
      t111 = t3 * s12
      t112 = t24 * s13
      t113 = t96 * s16
      t114 = ddreal(12)
      t115 = t24 * s25
      t116 = -t113 + t115
      t117 = t2 * s16
      t118 = ddreal(9) * s16
      t119 = (t117 - t84) * s12
      t120 = t79 * s16
      t121 = t96 * s25
      t122 = t95 * s26
      t123 = t96 * s26
      t124 = (s13 + s34 - s46) * s25
      t125 = -s13 + s46
      t126 = s16 * t125
      t127 = s16 - s34
      t128 = s16 * t127
      t129 = ddreal(10)
      t130 = t95 * s16
      t131 = t129 * s25
      t103 = t103 * s16
      t132 = s25 - t101
      t133 = s16 + s25
      t134 = s16 * t133
      t135 = t2 * s13
      t136 = t135 * s16
      t137 = t46 * s16
      t138 = t3 * (-t87 + t90)
      t139 = t2 * s26
      t140 = s12 * s16
      t141 = s16 + s25 - s34
      t142 = t129 * s26
      t143 = t3 * t98
      t144 = t5 * (-(-(-t137 + t115) * s26 - s12 * t132 - (t86 - t81) * 
     &s25 + t134 * t96 - t136) * s26 + ((t123 + t4) * s16 - t120 - t138)
     & * s34 + (-s16 * (s16 - s34 + t139 + t84) + t92) * s46 - (-(-t99 +
     & t4) * s56 + (-s12 - t122 - t121) * s25 + t126 * t2 + t24 * t141 *
     & s16 + t3 * (t124 + t140) + t142 * s16) * s56 + t105 - t106 - t107
     & - t108 * (s12 + t117 - t81) + t143)
      t145 = t10 * propW34
      t89 = -t8 * propW34 * (t3 * (-t24 * t98 + t3 * (-t110 * t116 - t10
     &5) + (s25 * (t113 - t112 + t111) + t104 + t103) * s16 + (s13 * t11
     &6 + (-s25 * t114 + t118) * s25 + t87 * t95 + t119) * s26 + ((-t122
     & - t121 - t101) * s16 + t24 * t90 + t120) * s34 + (-t115 * t58 + s
     &16 * (s16 - s34 + t123 + t115)) * s46 + ((t130 - t86) * s56 + s25 
     &* (-t131 + t130) + (ddreal(17) * s16 - ddreal(14) * s25) * s26 + t
     &124 * t24 + t126 * t96 + t128 * t46 + t119) * s56) + t144) + t145 
     &* (t3 * (-t110 * t85 + t98) + (-t4 * s16 + s25 * (s12 + t86) + t87
     & + t89) * s26 - (t3 * (t87 + t90) - (s25 + s13 - s26 - s12) * s16)
     & * s34 + (s16 * (s16 - t91 - t84) + t92) * s46 + (t102 * s56 + (t1
     &00 - t99) * s26 - t3 * (s16 + s34 - s46) * s25 + t96 * t97 + t89 +
     & t93 - t94) * s56 + t105 + t106 + t107 - t109)
      t105 = epinv2 * t7
      t116 = s25 - s13 + s16 - s34 + s46 + s56 + s12
      t119 = t3 * s26
      t144 = t116 + t119
      t146 = intHLs16s25s26s34s56x1121D2eps1()
      t147 = t3 * epinv
      t148 = -ddreal(1) + t147
      t149 = -s16 + s25
      t150 = t3 * t72
      t151 = t2 * s56
      t152 = t148 * s12
      t153 = intHLs160000x0111D0eps1()
      t154 = intHLs160000x0112D2eps1()
      t155 = intHLs160000x0121D2eps1()
      t156 = intHLs16s25s26s34s56x1112D4eps1()
      t157 = -t24 + t57
      t158 = intHL0s25s260s56x1011D2eps0()
      t159 = intHLs16s25s26s34s56x1221D4eps1()
      t160 = t24 * epinv
      t161 = -t2 + t160
      t162 = s16 + s26 - s34
      t163 = t46 * s56
      t164 = t3 * s56
      t165 = intHs16s25s26s34s56x1211D4eps1()
      t166 = s25 - s13 + s16 - s34 + s46
      t167 = t3 * t166
      t168 = t24 * s26
      t169 = s12 + t168 + t151 + t167
      t170 = intHL0s25s26s34s56x1110D2eps0()
      t171 = intHs160s26s34s56x1011D2eps0()
      t172 = intHs160000x0121D2eps1()
      t173 = t46 * s26
      t174 = s56 * t96 + t166 * t2 + s12 + t173
      t175 = intHL0s25s260s56x1010D0eps0()
      t176 = s12 + s13
      t177 = s25 - s13 + s16 - s34 + s12
      t178 = s25 - s34
      t179 = t2 * t178
      t180 = s45 + s12
      t181 = ddreal(9) * s26
      t182 = (s12)**2
      t183 = (s13)**2
      t184 = (t177 + t119) * s45
      t185 = t97 + t104
      t186 = (s56)**2
      t187 = t110 + t186
      t188 = t24 * t187 + t3 * t185 - s16 * t176 + (-t135 + t101) * s25 
     &- (s12 - t86 + t112 - t101) * s26 + (-t173 - t115 + t135 - t101) *
     & s34 - (-s16 - s45 - t168 - t179 + t81 - s46) * s46 + (t24 * t125 
     &+ t178 * t46 + t117 + t180 + t181) * s56 - t182 + t183 + t184
      t189 = intHL0s25s260s56x1020D2eps0()
      t190 = intHs16s25s26s34s56x1141D6eps1()
      t191 = -s12 + s16 + s26 - s34 + s56
      t192 = intHLs16s25s26s34s56x1111D2eps0()
      t193 = intHs160s26s34s56x1022D4eps1()
      t194 = t67 - t111
      t195 = intHs160s26s34s56x1031D4eps1()
      t196 = intHLs16s25s26s34s56x1112D2eps0()
      t197 = s25 - s16 + s34 + s56
      t198 = t197 + t111
      t199 = intHLs16s25s26s34s56x1113D4eps1()
      t200 = s25 + s56
      t201 = s46 + s45 - s12
      t202 = t3 * t201
      t203 = intHLs16s25s26s34s56x1212D4eps1()
      t204 = t2 * t201
      t205 = t96 * epinv
      t206 = intHLs16s25s26s34s56x1123D6eps0()
      t207 = t3 * t200
      t208 = s26 + t204 + t207
      t209 = intHLs16s25s26s34s56x1132D6eps0()
      t210 = t2 * t53 + t58
      t211 = intHLs16s25s26s34s56x1222D6eps0()
      t212 = s25 + s34 + s26
      t213 = t24 * t201
      t214 = t212 + t213 + t151
      t215 = intHLs16s25s26s34s56x1122D4eps1()
      t216 = ddreal(8)
      t201 = -s56 * t216 - t201 * t95 + s25
      t217 = intHLs16s25s26s34s56x1121D4eps1()
      t218 = intHL0s25s260s56x1011D2eps1()
      t219 = intHL0s25s26s34s56x1110D2eps1()
      t220 = intHs160s26s34s56x1011D2eps1()
      t221 = intHL0s25s260s56x1010D0eps1()
      t222 = intHL0s25s260s56x1020D2eps1()
      t223 = intHs16s25s26s34s56x1231D6eps0()
      t224 = s25 - s16 + s34 - s26 - s56
      t225 = t224 + t111
      t226 = intHs16s25s26s34s56x1222D6eps0()
      t227 = -s16 + s25 + s34
      t228 = s26 + s56
      t229 = t3 * t227
      t230 = t24 * s12
      t231 = intHLs16s25s26s34s56x1113D4eps0()
      t232 = intHLs16s25s26s34s56x1141D6eps0()
      t233 = intHLs16s25s26s34s56x1212D4eps0()
      t234 = intHLs16s25s26s34s56x1231D6eps0()
      t235 = intHLs16s25s26s34s56x1311D4eps0()
      t236 = intHLs16s25s26s34s56x1321D6eps0()
      t167 = t168 + t151 + t167
      t237 = t18 * s34
      t238 = t18 * s25
      t239 = t18 * s16
      t240 = t148 * s13
      t241 = s46 + s56
      t242 = t160 * s26
      t243 = -t148 * s12
      t244 = t141 * t167
      t245 = intHs16s25s26s34s56x1221D4eps1()
      t246 = s13 + t111
      t247 = t96 * s12
      t248 = t2 * t149
      t249 = s13 + s25
      t250 = -s12 + s16
      t251 = t3 * t250
      t252 = t2 * s12
      t253 = s13 + s25 - s46
      t254 = s16 - s34 - s12
      t255 = t24 * s34
      t256 = s13 - s46 - s56
      t257 = t3 * s34
      t258 = t2 * t133
      t259 = t216 * s26
      t260 = s16 * s26
      t261 = s16 * s56
      t262 = s16 + s56
      t263 = s25 * s26
      t264 = s16 + s46
      t262 = -(t259 + t100) * s34 + t2 * ((s46 + s25) * s25 + t261) - t2
     &4 * ((s46 - s34) * s34 + t182 - t260) + t3 * s16 * t264 + t46 * (-
     &s34 * t262 + t263) + t96 * s25 * t262 + s12 * (t24 * t256 + t257 -
     & t258 - t259) + s13 * (t255 - t4 - t101)
      t265 = s12 * t79
      t266 = (-t246 + s16) * s16
      t267 = -t2 * t265
      t268 = intHs16s25s26s34s56x1212D4eps1()
      t269 = s13 + s16
      t270 = s16 + s26
      t271 = -t2 * t270 + s12 - t115 + t81
      t272 = s16 - s26
      t273 = (s12 - t117 + t81) * s25
      t274 = (t269 + t4) * s26
      t275 = (-t3 * (s25 - s34 - s12) - t272) * s46
      t276 = t24 * (t110 - t97 + t186 - t104 - t265)
      t277 = t151 * (-s25 - s13 + s34 + s46 + t139 + t111)
      t278 = ddreal(9) * s12
      t279 = t278 * s26
      t280 = (t252 + t81) * s16
      t281 = intHs16s25s26s34s56x1211D2eps1()
      t282 = s13 - s16
      t283 = s12 + s25
      t284 = t3 * t282
      t285 = -t252 + t284
      t286 = s25 - s12
      t287 = -t286 + t119
      t288 = s13 - s16 + s34 - s46 - s12
      t289 = s12 - s13 + s16
      t290 = t289 * s12
      t291 = (t282 - s25) * s25
      t292 = t24 * t110
      t293 = -s34 + s46 + s56
      t294 = t293 * t283
      t295 = (-t282 + t111 + s25) * s25
      t296 = t119 * t283
      t297 = intHs16s25s26s34s56x1122D4eps1()
      t298 = t230 - t101
      t299 = -s13 + s25
      t300 = t2 * t299 - s12 + t113 + t259
      t301 = t3 * t270
      t302 = t2 * s34
      t303 = t230 + t302 - t301
      t304 = t96 * s34
      t305 = -t3 * (s25 - s13 + s46 - s12) - t173 - t99 + t304
      t306 = ddreal(11)
      t307 = t95 * s13
      t308 = -ddreal(15) * s25
      t309 = t216 * s13
      t310 = t46 * s12
      t311 = t216 * s34
      t312 = t95 * s12
      t313 = s16 * (-t176 + s16)
      t314 = t3 * t110
      t315 = t95 * t97
      t308 = t104 * t216 + t24 * t313 + t46 * t265 + (s16 * t306 - t252 
     &- t307) * s25 - (s12 * t114 + s13 - t118 + t308) * s26 + (-s16 * t
     &114 - ddreal(17) * s26 + t111 + t308 + t309) * s34 + (s26 - t311 +
     & t100 - t310 + t99) * s46 + (s25 * t306 - s34 * t114 + s56 + t119 
     &+ t137 - t312) * s56 + t314 + t315
      t316 = t3 * (-(-t2 * t250 - t299) * s26 + t186 + t313)
      t317 = t24 * (-t110 - t265)
      t318 = t2 * t104
      t319 = intHs16s25s26s34s56x1121D2eps1()
      t320 = t8 * t25 + t10
      t321 = -t299 - t139 - t101
      t322 = t3 * t127
      t323 = s25 - s13 + s46 + t139 + t322
      t324 = -s25 - s13 - s16 + t247
      t325 = s13 - s26
      t326 = t3 * (s25 + s16 - s12)
      t327 = -t325 + t326
      t328 = -s25 - s16 + s34 + s26 + t111
      t329 = s26 + s12
      t330 = t2 * t329
      t331 = -s25 - s13 - s16 + s34 + s46 + t330
      t332 = (-s13 + t251 + s25) * s25
      t333 = t3 * (t110 - t265)
      t334 = s13 * s16
      t335 = s12 * s13
      t336 = t250 * s25
      t337 = (-s25 + s13 + s12 - t117) * s26
      t338 = (-s12 + s16 + s26 - s34) * s46
      t339 = (-t323 - s56) * s56
      t340 = (-t321 - s34) * s34
      t341 = intHs16s25s26s34s56x1131D4eps1()
      t342 = t216 * s25
      t343 = t110 + t104
      t313 = t3 * t313
      t301 = -t2 * (-t97 - t186 - t265) + t24 * t343 + (t113 - t135 - t1
     &11) * s25 + (t342 + t137 - t247 - t81) * s26 + (s12 - t142 - t100 
     &- t137 + t112) * s34 + (t2 * t286 - t255 + t301) * s46 + (-t3 * (s
     &13 - s46 + s12) + t58 * t95 + t113 - t311) * s56 + t313
      t311 = intHs16s25s26s34s56x1211D2eps0()
      t344 = s13 - s16 + s34 - s46
      t345 = t3 * t228
      t346 = s16 - s34 + s46
      t347 = -t71
      t348 = s25 - s13 + s16 - s34 + s46 + s56
      t349 = s34 + s56
      t350 = epinv * s25
      t351 = t71 * s13
      t352 = epinv * s56
      t353 = epinv * s46
      t354 = t2 * (s56 - t353)
      t355 = t24 * (s26 - t352)
      t356 = t173 * epinv
      t357 = t24 * t178
      t358 = s25 * s46
      t359 = s34 * s46
      t360 = t241 * s56
      t361 = s56 * t178
      t362 = s26 * s56
      t363 = t347 * s12
      t364 = (-t348 - t119) * t349
      t365 = intHs16s25s26s34s56x1311D4eps1()
      t366 = (t166 + t139 + t164 + s12) * s12
      t167 = t228 * t167 + t366
      t367 = intHs16s25s26s34s56x1113D4eps0()
      t368 = (s16 + t84) * s26
      t369 = intHs16s25s26s34s56x1131D4eps0()
      t370 = intHs16s25s26s34s56x1212D4eps0()
      t371 = intHs16s25s26s34s56x1120D2eps1()
      t372 = -s12 - s13 + s25 + s45 + s46 + s56
      t373 = intHs16s25s26s34s56x1311D4eps0()
      t374 = intHs160s26s34s56x1012D2eps0()
      t375 = intHs160s26s34s56x1013D4eps0()
      t376 = intHs160s26s34s56x1021D2eps0()
      t377 = intHs160s26s34s56x1022D4eps0()
      t378 = intHs160s26s34s56x1031D4eps0()
      t379 = intHs160s26s34s56x1021D2eps1()
      t380 = t18 * s12
      t381 = s12 + s34 + s56
      t382 = s12 + s34
      t383 = t2 * t382
      t384 = intHs16s25s26s34s56x1132D6eps0()
      t385 = s16 + s26 + s56
      t386 = t3 * t385
      t387 = t383 - t386
      t388 = intHL0s25s260s56x1022D4eps0()
      t389 = t3 * t60 + t228
      t390 = intHs16s25s26s34s56x1220D4eps0()
      t391 = -t3 * (s13 - s46 - s45 + s12) - s16 + s25 - s26 + s34 + s56
      t392 = intHs160000x0111D0eps0()
      t393 = t166 + t345
      t394 = intHs160000x0112D2eps0()
      t395 = intHs16s25s26s34s56x1111D2eps0()
      t396 = intHLs16s25s26s34s56x1123D6eps1()
      t397 = intHs160000x0211D2eps0()
      t398 = intHs16s25s26s34s56x1211D4eps0()
      t399 = intHLs16s25s26s34s56x1132D6eps1()
      t400 = intHLs16s25s26s34s56x1222D6eps1()
      t401 = intHs160000x0121D2eps0()
      t402 = intHs16s25s26s34s56x1121D4eps1()
      t403 = -t2 + epinv
      t404 = intHs16s25s26s34s56x1221D4eps0()
      t405 = intHs16s25s26s34s56x1121D2eps0()
      t406 = intHs16s25s26s34s56x1122D4eps0()
      t407 = intHs160s26s34s56x1012D2eps1()
      t408 = intHs160s26s34s56x1013D4eps1()
      t409 = intHs16s25s26s34s56x1210D2eps1()
      t410 = intHs16s25s26s34s56x1310D4eps1()
      t411 = intHL0s25s260s56x1031D4eps1()
      t412 = -s12 + s25 + s26 + s45 + s46 + s56
      t413 = intHL0s25s260s56x1022D4eps1()
      t414 = intHs16s25s26s34s56x1220D4eps1()
      t415 = intHs160000x0111D0eps1()
      t416 = intHs160000x0112D2eps1()
      t417 = intHs160000x0211D2eps1()
      t418 = intHs16s25s26s34s56x1312D6eps1()
      t419 = intHs16s25s26s34s56x1321D6eps1()
      t420 = intHLs16s25s26s34s56x1312D6eps1()
      t421 = intHLs16s25s26s34s56x1114D6eps1()
      t422 = intHs16s25s26s34s56x1312D6eps0()
      t423 = intHs16s25s26s34s56x1321D6eps0()
      t424 = intHLs16s25s26s34s56x1312D6eps0()
      t425 = intHLs16s25s26s34s56x1114D6eps0()
      t426 = intHLs16s25s26s34s56x1213D6eps0()
      t427 = t227 + t202
      t428 = intHLs16s25s26s34s56x1211D2eps1()
      t429 = s34 - s46 - s45 + s12 - t164
      t430 = intHLs16s25s26s34s56x1311D4eps1()
      t431 = intHLs16s25s26s34s56x1122D4eps0()
      t432 = intHLs16s25s26s34s56x1221D4eps0()
      t433 = intHLs16s25s26s34s56x1131D4eps0()
      t434 = intHLs16s25s26s34s56x1121D2eps0()
      t435 = intHLs16s25s26s34s56x1211D2eps0()
      t436 = intHLs16s25s26s34s56x1213D6eps1()
      t437 = intHs16s25s26s34s56x1231D6eps1()
      t438 = intHs16s25s26s34s56x1222D6eps1()
      t439 = intHs16s25s26s34s56x1132D6eps1()
      t440 = intHs16s25s26s34s56x1130D4eps0()
      t441 = intHs16s25s26s34s56x1112D4eps0()
      t442 = t3 * t329
      t166 = t166 + t442
      t443 = intHs16s25s26s34s56x1112D4eps1()
      t444 = t24 * t116
      t445 = intHs16s25s26s34s56x1121D4eps0()
      t446 = s56 - s12
      t447 = intHs16s25s26s34s56x1120D2eps0()
      t448 = intHs16s25s26s34s56x1210D2eps0()
      t449 = intHs16s25s26s34s56x1310D4eps0()
      t450 = intHL0s25s260s56x1031D4eps0()
      t451 = intHL0s25s260s56x1012D2eps1()
      t452 = intHL0s25s260s56x1013D4eps1()
      t453 = intHs160s26s34s56x1020D2eps1()
      t454 = intHs16s25s26s34s56x1130D4eps1()
      t455 = intHs16s25s26s34s56x1141D6eps0()
      t456 = ddreal(1) / s25
      t457 = t190 * t18
      t458 = t457 - t455
      t459 = t421 * t18
      t460 = t459 - t425
      t461 = t52 * t18
      t462 = t461 - t232
      t463 = -t36 + t29
      t464 = ddreal(1) / (gw)**2
      t465 = ddreal(1) / (ecossin)**2
      t466 = t458 * t191
      t467 = t7 * t464
      t468 = t467 * t456 * t465 * propW34
      t469 = t54 * t18
      t470 = t468 * t12
      t471 = t470 * t53
      t472 = t439 * t18
      t473 = t437 * t18
      t474 = -t473 + t223
      t475 = t420 * t18
      t476 = t475 - t424
      t477 = t33 * t18
      t478 = t55 * t18
      t479 = (-t477 + t30) * (-t270 + t383)
      t480 = (-t472 + t384) * t387
      t481 = ddreal(1) / t144
      t482 = ddreal(1) / t67
      t483 = ddreal(1) / t27
      t484 = t418 + t419
      t485 = t191 * t378
      t313 = t367 * t320 * (-t24 * (t97 + t104 + t368) + t267 - t313 + (
     &s12 - t137 + t112) * s25 + (t216 * t58 - s12 - t112 + t137) * s34 
     &+ (-t357 + t252 - t101) * s46 + (-t178 * t96 - t272 * t3 + s56 + t
     &125 + t230) * s56 + t173 * s12) * t481
      t167 = t468 * (t12 * (t18 * (t208 * t396 + t210 * t399 + t427 * t4
     &36) + t483 * (-t13 * t58 + t412 * t450 + t60 * t61) + t48) + t26 *
     & (t18 * (-t198 * t34 - t283 * t484) + t481 * (t169 * t397 + t28 * 
     &(epinv * (t216 * (-s34 * t58 + t263) + t24 * (t97 + t104 - t359 + 
     &t358) + t3 * ((t241 + t179 + t119 + s16) * s16 - t362) + t96 * t36
     &1 + s13 * (s56 - t357 - t101) - t360) + (t3 * (s25 + s46 + t239) +
     & t351 - t350 - t237 + t354 + t355 - t356 + t363) * s12 - t364) + t
     &301 * t369 + t341 * (t3 * ((-t323 - s56) * s56 - t104 + t182 + t33
     &4 - t335 - t336 + t337 - t338 - t87) - t292 - t3 * t321 * s34 + ep
     &inv * t301) + t365 * (-epinv * t167 + t290 + t294 + t295 + t296) -
     & t167 * t373 + t393 * t394) + t482 * (s12 * t375 - t372 * t449 + t
     &440 * t66 - t485)) - t313)
      t301 = t191 * t195
      t313 = epinv * t481
      t74 = t468 * (t26 * (t313 * (t169 * t417 + t393 * t416) + (-s12 * 
     &t408 + t372 * t410 - t454 * t66 + t301) * t482 * t18) + (t17 * t58
     & - t411 * t412 - t452 * t60) * t18 * t483 * t12) + t167 + t468 * (
     &t12 * (t154 * epinv + t53 * (t478 + t235 - t236) - t70 * (epinv * 
     &(s25 - t74 - t73) + t76 + t77) - t199 * (-t160 * t53 + t200 + t202
     &) - t206 * t208 - t209 * t210 + t476 * (s34 + s46 + s45 - s12) - t
     &426 * t427 - t430 * (-epinv * t53 - s12 + s34 + s45 + s46) - t433 
     &* (s25 - t74 - t204)) + t26 * (t198 * t31 + t225 * t474 + t283 * (
     &t422 + t423) + t479 + t480))
      t167 = ddreal(1) / s16
      t198 = (t483)**2
      t208 = t158 * t483
      t210 = (-t221 + t222) * t18
      t323 = t210 * t198
      t412 = t323 * t481
      t427 = t438 * t18
      t486 = t456 * t464
      t487 = t486 * t465
      t488 = t392 - t395
      t489 = -t447 - t448
      t490 = t194 * t377
      t491 = t379 * (s16 - s34 + s26 + s56 - t380)
      t492 = t400 * t18
      t493 = -t388 * t389 * t483 - t492 * t214
      t346 = t481 * (-t311 * ((t122 + t151) * s56 + t3 * (s26 * t346 + s
     &56 * t346) + t292 + (-t344 + t139 + t164 + s12) * s12 + s13 * (s25
     & - t345) + (-s16 + s34 - s46 + s56 - s25) * s25) + t370 * (-t3 * (
     &t87 - t275 - t273 + t274) + t276 + t277 + t279 - t3 * t271 * s34 +
     & t280)) * t320
      t256 = t468 * (t12 * t493 + t26 * (t481 * (t169 * t398 - t174 * t4
     &01 + t23 * (-t182 * t3 + s12 * (t256 * t3 - t141 - t168) + t244) -
     & t262 * t404 - t268 * (t119 * (t149 + t252) + epinv * (t3 * (s34 *
     & t271 - t273 + t274 - t275 + t87) - t276 - t277 - t279 - t280) + (
     &-s13 + t230 + s25) * s25 + (t3 * (s16 + s26 - s12) - s13 - s34) * 
     &s34 + t241 * (t227 + t252) - t266 + t267) - t406 * t308 + t488 * t
     &393 + t402 * (-s12 * t403 + t2 * t348 + t173 + t352) + t405 * (-s2
     &6 * t324 + (-t327 + s34) * s34 - s46 * t328 + (-t331 - s56) * s56 
     &+ t266 + t332 - t333) - t443 * (epinv * t166 - t259 - t444) + t445
     & * t446) + t482 * ((t374 + t376) * s12 - t371 * (epinv * t372 + s1
     &2 + s13 + s16 + s26 - s34 - s45 - s46) + t372 * t489 + t390 * t391
     & + t490 + t491)) + t346)
      t271 = t18 * t453
      t274 = s12 * t407
      t275 = t68 * t69
      t276 = (t271 - t65) * t482
      t277 = t18 * t451
      t279 = t19 * t18
      t280 = t69 * t218
      t346 = t53 * t56
      t57 = t483 * (t58 * (-t279 + t15 + t21) + t60 * (-t277 + t59 + t62
     &) - t63 * (s25 + s46 + s26 + s56 + s45 - s12 - t64) - t280) - t192
     & + t47 - t49 + t51 + t346 * (-ddreal(1) + t57)
      t324 = t319 * t320 * (epinv * (s26 * t324 - (-t327 + s34) * s34 + 
     &s46 * t328 - (-t331 - s56) * s56 - t266 - t332 + t333) - t182 + t8
     &7 - t334 + t335 + t336 - t337 + t338 - t339 - t340 + t314) * t481
      t57 = t468 * (t12 * t57 + t26 * (t481 * (t35 * ((t3 * (-epinv * t2
     &41 + s26) + s46 + s56 - t237 + t238 + t239 + t240 - t242 + t243) *
     & s12 + t244 * epinv) - t245 * (t3 * (-t110 + t97) + epinv * t262 +
     & (s16 + t247 - t81) * s25 + (s13 + t247 + t248) * s26 - (t249 - t1
     &39 - t251 + s34) * s34 + (-t162 + t252 + t84) * s46 - (t254 * t3 +
     & s56 + t139 - t253) * s56 - t266 + t267) + t281 * (epinv * (-t186 
     &* t2 + s26 * t285 + s34 * t287 - t287 * s46 - (-t288 * t3 + s25 + 
     &t122) * s56 - t290 - t291 - t292) + t290 + t294 + t295 + t296) - t
     &297 * (epinv * t308 + s25 * t298 + s34 * t300 + s46 * t303 + s56 *
     & t305 - t316 + t317 - t318) + epinv * t393 * t415) + t482 * (t18 *
     & (t372 * t409 - t274) + t275 + t276 * t66)) - t324)
      t66 = t175 - t189
      t194 = t193 * t194
      t244 = t413 * t18
      t174 = t468 * (t12 * (t483 * (-t169 * t170 * t481 + t244 * t389) +
     & t481 * t188 * t66 * t198) + t26 * (t482 * (-t18 * (t391 * t414 + 
     &t194) + t169 * t171 * t481) + t313 * (t165 * t169 - t172 * t174)))
      t57 = t174 + t57 + t256 + t487 * (t105 * t89 * t167 * t483 * t482 
     &* t481 + (t12 * ((t153 - t155) * epinv - t146 * (-epinv * (-t149 +
     & t151 + t150) + s45 + s46 + s56 + t152) - t18 * t217 + t156 * t157
     & + t159 * (-epinv * (t24 * t72 + t162 + t163) + s12 * t161 + s34 +
     & t73 + t164) + t203 * (-t205 * t53 + t197 + t204) + t211 * t214 + 
     &t215 * (epinv * t201 + t24 * t53 + s25 + s26) + t431 * t201 - t428
     & * (epinv * t429 - s12 + s34 + s45 + s46) - t435 * t429 + t432 * (
     &-t162 - t163 - t213) - t434 * (t149 - t151 - t202) - t208 + t412 *
     & t188) + t26 * (t78 * t482 + (t427 - t226) * (-t228 + t230 + t229)
     & - t441 * t166 * t481) + t481 * (-t12 * t219 * t483 + t26 * t220 *
     & t482) * t169 * t69) * propW34 * t7)
      t73 = t196 * t53
      t89 = ddreal(8) / ddreal(9)
      t163 = -ddreal(5) / ddreal(9)
      t166 = ddreal(4) / ddreal(9)
      t169 = ddreal(2) / ddreal(3)
      t174 = ddreal(2) / ddreal(9)
      t188 = ddreal(1) / ddreal(3)
      t197 = ddreal(1) / ddreal(9)
      t201 = epinv2 * t167
      t42 = t483 * (-t1 - t14 + t42 - t201)
      t202 = (t18 * (t268 + t297) - t23 - t370 - t406 - t44) * t26
      t204 = t199 * epinv
      t213 = t483 * t40
      t214 = (t18 * (t438 + t439 + t418) - t226 - t28 - t384 - t422) * t
     &26
      t256 = t33 + t34
      t262 = t256 * t18
      t266 = t24 * t26 * (-t262 + t30 + t31)
      t308 = -t46 * t38
      t324 = (-t146 - t428) * epinv
      t327 = -t459 + t425
      t328 = t327 * t12
      t331 = t46 * (-t38 + t328)
      t332 = s25 + s26 - s34 + s56
      t333 = -s13 + s25 + s26 - s34 + s56
      t348 = s26 - s34 + s56
      t352 = s56 * t58
      t372 = -t2 * t27 * s34 + t3 * (t104 + t352) - s13 * (s25 + s16 + s
     &26 + s56 - t257) + t348 * s16 + (-s46 + s26 - s12) * s25 + t110 + 
     &t186 + t97
      t389 = s16 - s26 - s56
      t391 = s46 + s56 + s12
      t393 = s13 * t389 + t332 * (t391 + t119)
      t429 = t3 * t178
      t493 = s16 + t429
      t494 = s34 - s56
      t495 = epinv * t493
      t496 = s16 + s25 - s34 + s56
      t497 = s25 - s34 + s26
      t498 = s26 * s34
      t499 = s25 * s34
      t497 = -t3 * (-s56 * t497 + t499) - s13 * t496 + s16 * (t497 + t16
     &4) + s56 * t391 + t104 + t263 + t97 - t498
      t500 = s25 - s13 - s34 + s56
      t501 = s16 + s46 + s12
      t502 = t3 * t500
      t503 = t501 + t139 + t502
      t504 = epinv * t503 - s12 + s13 - s16 - s25 + s34 - s46 - s56 - t1
     &19
      t505 = s25 - s13 + s16 + s12
      t506 = t24 * t325
      t507 = t58 * s46
      t255 = -t104 * t24 - t2 * t97 + (-t250 + t135) * s25 + (-t250 - t8
     &4 + s26) * s26 + (t100 - t506 + t101) * s34 + (-s16 + s46 + s26 + 
     &s12 + t255 - t4) * s56 + t507 + t81 * s16
      t508 = (-s46 - s56 - t505 - t119 + s34) * s34
      t509 = s26 + s34 + s56
      t510 = t3 * t509
      t511 = -s46 - s12
      t512 = s16 * s34
      t513 = t27 * t144
      t514 = t2 * t500 + t168 + t501
      t515 = -t3 * (-(s46 + s26 + s12) * s25 - s34 * t332) + s13 * (s16 
     &- t257) + (s25 + s34 - s26 - s56) * s16
      t516 = s26 - s34
      t517 = s13 * t27
      t518 = t3 * t352 + (s16 - s34 + s46 + s12) * s25 + s26 * t516 - s3
     &4 * s56 + t186 + t97 - t517 + t4 * s26
      t519 = ddreal(1) - t9
      t520 = -t2 + t147
      t521 = s12 + s16 + s26 - s34 - s45 - s46 + s56
      t522 = s26 - s56
      t523 = s25 + s26 + s34 - s56
      t524 = intHLs160000x0122D4eps0()
      t525 = t137 * propW16
      t526 = t525 * t519
      t527 = t8 * (-t3 + t5 - t526) + t10
      t528 = t526 * t8
      t529 = intHLs160000x0111D2eps0()
      t530 = t3 * t299
      t531 = s25 + s46 + s12
      t532 = t2 * t110
      t533 = t80 * s16
      t534 = s13 * t132 + (t113 - t84) * s25 + t533
      t535 = -s25 + t117
      t536 = s13 - s46 - s12
      t537 = t535 * s26
      t538 = (s13 + t101) * t97
      t94 = (-s34 * t132 - t94) * t58
      t539 = t96 * s13
      t540 = ddreal(13)
      t541 = t540 * s16
      t542 = s25 + s13 - s12
      t543 = t24 * t97
      t544 = t2 * t87
      t545 = t216 * s16
      t535 = (-t535 * s26 - t534) * s26
      t546 = t513 * t101 * propW16 * t67
      t547 = t513 * t525 * t67
      t548 = intHLs160000x0122D4eps1()
      t549 = t1 + t14
      t550 = (t397 + t401) * t514
      t551 = t320 * t481
      t552 = -t175 + t189
      t553 = t69 * t219
      t554 = t332 * t483
      t555 = epinv * t153
      t556 = t156 * t520
      t557 = t483 * t481
      t558 = t557 * t333
      t559 = t220 * t69
      t560 = t228 * t281
      t561 = t18 * t409
      t562 = t224 * t320
      t94 = -t562 * t390 - t201 * (t8 * (t24 * (t98 + t535 + t94 - t109 
     &- t538) - t3 * ((-s25 - s13 + s46 + s12 + t137 + t139 - t257 + s56
     &) * s56 + (-t539 + t252) * s16 + (t118 + t81) * s25 + (-t79 + t541
     & - t4) * s26 + (t542 - t137 - t139 + s34) * s34 + (t516 + t117) * 
     &s46 + t544 - t543 + t314) * s56 + t547 - t5 * (((t536 - t99 - t139
     & + t257 - s56) * s56 + t3 * (-t110 + t97 - t87) + (t135 - t111) * 
     &s16 - (s13 + t113) * s25 + (t542 - t545) * s26 - (s34 + t79 - t99 
     &- t139) * s34 + (-t516 - t101) * s46) * s56 - t109 - t538 + t94 + 
     &t98 + t535 + t546)) + t10 * ((-t534 - t537) * s26 + ((-t536 + t139
     & + t429 + s56) * s56 + t3 * (t110 + t97) - (s13 + t117) * s25 - (t
     &79 - t4 + t101) * s26 + (t79 - t139 - t84 + s34) * s34 + t516 * s4
     &6 + t334) * s56 - t109 + t98 - t538 + t94)) * t481 * t483
      t109 = t467 * t465 * propW34
      t94 = t109 * (t26 * (t456 * (t481 * (-t165 * t504 - t35 * (epinv *
     & t497 + t508) + t393 * t405 - t559 * t503 * t482 + t560 * epinv * 
     &t333) + t482 * (t65 - t78)) - t561 * t482) + t456 * (t482 * t94 + 
     &t551 * t297 * (-epinv * t255 + t508)) + t456 * (-t203 * (epinv * t
     &521 + s25 + s26) - t211 * t523 - t215 * (s25 + s26 - t64) - t233 *
     & t521 + t431 * t60 + t434 * t58 + t435 * t494 - t73 - t555 - t556 
     &+ t558 * (t554 * t552 + t553)) * t12)
      t118 = t172 + t417
      t508 = t221 - t222
      t521 = t311 * t228
      t534 = t171 * t482
      t535 = t534 * t503
      t538 = t118 * epinv
      t542 = t18 * s56
      t563 = t48 * t527
      t527 = t154 * (epinv * t527 + t528)
      t564 = t18 * t414
      t393 = t109 * (t26 * (t482 * (-t275 * t456 + t447 + t448) + t456 *
     & (epinv * (-t245 * t515 + t319 * t393) + t333 * (-(t415 + t443) * 
     &epinv + t521) - t535 + t538 * t514) * t481 + t456 * t482 * (epinv 
     &* t379 + t18 * (t193 - t407)) * s34) + t456 * (t564 * t562 * t482 
     &+ t527 + t563) + t456 * (t483 * (t542 * (t21 + t22) + (t554 * t18 
     &* t508 + t170) * t481 * t333) + t192 - t47) * t12)
      t36 = t36 - t29
      t562 = epinv * t155
      t565 = t18 * t548
      t566 = t69 * intHLs160000x0111D2eps1()
      t567 = t26 * t36
      t568 = t456 * t465 * propW34 * t7
      t459 = t568 * (t464 * (t567 * s56 + (-t459 + t425) * (s25 - s56) *
     & t12) + ((t565 + t49 - t524 + t562) * s16 - t529 - t566) * t519 * 
     &propW16)
      t569 = t396 * t18
      t570 = -t569 + t206
      t571 = t399 * t18
      t572 = -t571 + t209
      t573 = -t475 + t424
      t574 = t436 * t18
      t575 = -t574 + t426
      t576 = t410 * t482
      t577 = t367 * t493
      t578 = t456 * t26
      t64 = t109 * (t456 * (t12 * t50 - t26 * t440) + t26 * (-t237 * t40
     &8 * t456 + t449) * t482 + t542 * t456 * t12 * (t17 + t19 + t20) * 
     &t483 + t578 * ((-t333 * t416 - t341 * t372 + t365 * t518) * epinv 
     &- t577 * t333) * t481) + t109 * (t26 * (t18 * (t456 * (s56 * t256 
     &+ t454) - t576) + t456 * (-(t30 + t31) * s56 + s34 * t375 * t482) 
     &+ (t28 * (t144 * t494 - t495 * t333) - t333 * t394 - t369 * t372 +
     & t373 * t518) * t456 * t481) + t456 * (t199 * (s25 - s56 - t64) - 
     &t231 * t60 + t494 * t573 + t58 * t572 + t575 * (t227 - t164) + t57
     &0 * (t522 + t84) + (-t13 - t15 - t16) * t483 * s56) * t12)
      t256 = -t18 * intHLs160000x0113D4eps1() + intHLs160000x0113D4eps0(
     &)
      t372 = ddreal(4) / ddreal(3)
      t518 = t372 * t568 * propW16 * s16 * t519 * t256
      t64 = t169 * t459 + t174 * t64 - t197 * (t393 + t94 + t468 * (t12 
     &* ((t146 * t58 + t428 * t494 - t346) * epinv + t492 * t523) + t26 
     &* (t481 * (-t23 * t497 - t268 * (epinv * (t2 * (t110 + t499) - t3 
     &* (-t186 + t104 - t263) + s13 * (t133 - t510) + s25 * t391 - t511 
     &* s26 - s56 * t511 - t97 + t512 + t123 * s56) + t513) + t333 * (-t
     &392 + t395 - t441) - t402 * t504 - t404 * t515 + t503 * (-t398 - t
     &445) + t550) + t482 * ((t374 + t376 - t377) * s34 + t371 * (s16 - 
     &s34 + s26 + s56 + t350) - t271)) - t12 * t549 * t483 * s56 - t551 
     &* (t255 * t406 + t370 * (-t3 * (-t186 + t104) - (-t176 + s25) * s2
     &5 + (s12 + t530) * s26 + (s16 + t4 - t81) * s34 + (t531 + t123 - t
     &81) * s56 + t334 + t507 + t532)))) - t518
      t94 = intHLs160000x0211D2eps0()
      t255 = intHLs16s25s26s34s56x1211D4eps0()
      t346 = s13 + s25 + s26 + s56
      t391 = t2 * t249
      t393 = t228 + t391
      t459 = t228 + t357 + t101
      t497 = t385 + t429
      t504 = t496 + t119
      t507 = epinv * s26
      t511 = t346 + t257
      t515 = s25 + s26 + s34 + s56
      t542 = -s46 + s26 + s56 - t307 + t84
      t579 = t3 * t346 + s34
      t580 = t391 + t345
      t581 = t3 + epinv
      t347 = t347 * s13
      t582 = s25 - s16 + s26 + s56 - t81
      t583 = -t161
      t584 = t583 * s13
      t585 = t24 * t249
      t510 = -s16 + t585 + t510
      t586 = intHLs160000x0211D2eps1()
      t587 = intHLs16s25s26s34s56x1211D4eps1()
      t588 = s13 - s34
      t589 = t3 * t249
      t590 = ddreal(1) + epinv
      t591 = t18 * s13
      t592 = s34 * t590 + t591
      t593 = s25 - s46
      t594 = t593 - t135
      t595 = s25 - s16 + s34 - s46 - t539
      t596 = -t3 * t148 * s13
      t597 = t389 + t429
      t598 = t249 * t18
      t599 = t228 + t589
      t600 = s25 + s34
      t601 = s25 * t27
      t602 = (t59 + t62) * t249
      t603 = t63 * (s26 + s56 + t598)
      t604 = (t244 - t388) * t599
      t605 = t428 * t592
      t521 = t468 * (t12 * (-t146 * (epinv * t582 + s13 + s25 + s26 + s5
     &6) + t159 * (t515 * t581 + t584) - t170 * t483 - t196 * t594 + t20
     &3 * (epinv * t595 + t135 + t229) + t211 * t510 + t215 * (epinv * t
     &542 + t345 + t585) + t233 * t595 + t431 * t542 - t434 * t582 - t56
     & * (epinv * t594 + s13 + s25) - t605) + t26 * ((t141 * (t35 + t319
     &) - t245 * t493 + t402 + t415 + t443) * epinv + t141 * t405 - t370
     & * t597 - t406 * t459 - t521 - t371 * (s25 + s16 - s34 + s56 - t50
     &7 + t119) * t482))
      t542 = t598 * t451
      t582 = t435 * t588
      t585 = t404 * t493
      t6 = t468 * (t12 * (t157 * (t156 + t217 + t587) + t432 * (t515 - t
     &112) + t483 * (-t280 - t158 - t542) + t47 + t555 + t582) + t26 * (
     &epinv * t165 + t141 * t23 + t482 * (t275 + t78 + t171) + t392 - t3
     &95 + t398 + t441 + t445 - t585)) + t521 + t468 * (t26 * ((-t268 * 
     &t597 - t297 * t459 - t560) * epinv + t482 * ((-t561 + t447 + t448)
     & * s26 + (t564 - t390) * (t496 + t139) + t559 + t276 * t504)) + t4
     &83 * (t12 * (t604 - t553 + t602 - t603) + t201 * (t8 * ((s16 * t6 
     &+ (t46 * propZ25 * t228 - t302 * propZ25 - t24 + t5) * s25 - t228 
     &* t46 + t257) * s16 + t601 * t25) + t10 * ((-t600 + s16) * s16 + t
     &601)) * t482) - t492 * t12 * t510 + (-s13 - s34 + s46 + s26 + s45 
     &+ t207) * t12 * (t210 + t175 - t189) * t198)
      t25 = (-t249 * t61 - t346 * t450) * t483
      t207 = t228 * t365
      t210 = t341 * t497
      t275 = t228 * t373
      t276 = t249 * t452 + t346 * t411
      t459 = t234 * t579
      t496 = t235 * t588
      t497 = t369 * t497
      t510 = t462 * t346
      t521 = t460 * t249
      t560 = s34 * (t18 * intHLs16s25s26s34s56x1411D6eps1() - intHLs16s2
     &5s26s34s56x1411D6eps0())
      t6 = t169 * t470 * (t560 + t510 + t521) - t174 * (t468 * (t12 * (-
     &(t154 + t155 + t586) * epinv + t206 * t393 + t209 * t580 + t236 * 
     &t511 + t70 * (t27 * t581 + t347) - t48 - t49 - t94 + t459 - t496) 
     &+ t26 * (-(t172 + t416 + t417) * epinv - t394 - t397 - t401 - t577
     & - t497) + t12 * t276 * t483 * t18) + t468 * (t12 * (-t18 * (t393 
     &* t396 + t399 * t580 + t511 * t55 + t54 * t579) - t199 * (-s25 * t
     &581 + t353 - t596) + t231 * (t593 - t112) + t430 * (s34 * t581 + t
     &591) + t433 * (t27 - t135) + t573 * (s25 + s13 - s16 + t257) + t57
     &5 * (-t127 + t589) + t25) + t26 * ((-t28 * t493 + t207 - t210) * e
     &pinv + t482 * ((t18 * t410 - t449) * s26 + t504 * (t18 * t454 - t4
     &40)) + t275))) - t188 * t470 * (t50 + t51 + t255) - t197 * t6
      t493 = t18 * (t203 + t215)
      t39 = t39 * t483
      t504 = (t13 + t15 + t16) * t483
      t577 = (t396 + t436) * t18
      t588 = -ddreal(1) + epinv
      t593 = t588 * t418
      t594 = t56 * epinv
      t11 = t11 * t196
      t595 = t34 * t18
      t597 = epinv * (t8 * (t203 + t215) + t10 * (t28 + t268 + t297 + t3
     &99 + t400 - t418 + t420 - t438 - t439)) + t12 * t199 + t8 * (t431 
     &+ t233 - t203 - t215 + t594) + t10 * (t438 + t439 + t418 - t420 - 
     &t422 + t424 - t399 - t400 + t406 - t384 + t367 + t370 - t226 + t20
     &9 + t211 - t28)
      t598 = t588 * t32
      t601 = t9 * t8
      t606 = ddreal(1) - epinv
      t607 = t588 * t396
      t608 = t588 * t436
      t609 = t588 * t33
      t610 = t588 * t34
      t611 = t588 * t438
      t612 = t588 * t439
      t613 = t8 * (-epinv * (t268 + t297) + t28 * t606 + t206 + t226 - t
     &367 - t370 + t384 - t406 + t422 + t426 + t593 + t607 + t608 + t611
     & + t612)
      t597 = -t114 * t8 * (t421 * t606 - t425 + t9 * (t610 + t30 + t31 +
     & t206 + t426 + t607 + t608 + t609)) + t2 * t601 * (t203 * t588 + t
     &215 * t588 + t233 - t268 - t297 + t431 + t594) + t216 * t613 - t24
     & * (t8 * (-t424 - t297 - t268 - t209 - t211 + t492 + t475 + t571) 
     &+ t10 * (t426 + t206 - t30 - t31 - t574 - t569 + t595 + t477)) - t
     &3 * t597 - t46 * (t9 * epinv * (t7 * (-t32 + t421) + t8 * (-t28 - 
     &t268 - t297 + t399 + t400 + t420 + t438 + t439)) + t9 * (t7 * (-t2
     &9 + t32 - t421 + t425) + t8 * (t593 + t28 + t209 + t211 + t226 - t
     &367 - t370 + t384 - t399 - t400 - t406 - t420 + t422 + t424 - t438
     & - t439))) - ddreal(18) * t601 * (t421 * t588 + t29 + t425 + t598)
     & + ddreal(24) * t8 * (t598 + t29) + ddreal(16) * t8 * (t610 + t30 
     &+ t31 + t609) + t10 * (-t18 * t203 - t18 * t215 + t196 + t233 + t2
     &68 + t297 + t431 + t594) + t11
      t598 = t3 * (t12 * (t18 * (-t436 - t420 - t399 - t400 - t396 + t39
     &) + t206 + t209 + t211 + t424 + t426 - t504) + t214) - t266 - t308
     & - t12 * (-t43 * t483 + t233 + t431 - t493) - t202
      t607 = t2 + epinv
      t613 = t13 * t483
      t614 = t428 * t590
      t615 = t211 - t492
      t616 = t12 * t615
      t617 = t2 * t616
      t19 = t24 * (t12 * (t18 * (t54 + t399) - t209 - t234) + t26 * (-t3
     &0 + t477)) + t3 * (t12 * (t18 * (t436 + t420 + t396) + t434 - t613
     &) + t26 * (t18 * (t439 + t34) - t31 - t384) + t12 * (-t215 * t590 
     &+ t483 * (t18 * (-t452 - t413 - t411 + t17) + t388 + t450 + t61) -
     & t581 * t70 - t206 - t236 - t424 - t426 - t430 - t431 - t433 + t47
     &8) + t26 * (-t18 * t28 + t367 + t406)) + t46 * (t12 * (t461 - t232
     &) + t38) + t12 * (t146 * (ddreal(1) + t147) + t483 * (t18 * (t451 
     &+ t63 + t19) - t15 - t21 - t59 - t62) + t607 * (-t159 - t203) - t1
     &96 - t233 - t432 + t435 - t594 + t614) + t26 * (t148 * t297 - t226
     & - t268 + t427) - t617
      t457 = t457 - t455
      t461 = t26 * t457
      t477 = t454 * t482
      t618 = t440 * t482
      t619 = (-t203 + t215 - t56) * epinv
      t620 = t319 * t590
      t621 = (t390 + t371 - t564) * t482
      t622 = t26 * (t226 - t427)
      t623 = t46 * t461
      t624 = -t2 * t622 + t24 * t26 * (t18 * (t437 + t439) - t223 - t341
     & - t384) + t3 * (t12 * (-t18 * (t436 + t420) + t204 + t231 + t424 
     &+ t426) + t26 * (t18 * (t418 + t419 + t33 + t34 - t477) - t245 - t
     &268 - t297 - t30 - t31 - t422 - t423 + t618)) + t12 * (t431 - t233
     & - t196 + t211 + t619 - t492) + t26 * (t405 + t621 + t620) + t623
      t625 = t588 * t419
      t626 = t18 * t419
      t611 = t611 + t226
      t612 = t437 * t588 + t223 + t341 + t384 + t612
      t627 = t8 * (t438 * t606 + t9 * t612 - t226)
      t628 = t190 * t588 + t455
      t629 = t8 * (t625 + t30 + t31 + t245 + t268 + t297 + t422 + t423 +
     & t593 + t609 + t610)
      t11 = -t114 * t627 + t2 * t9 * (t611 * t7 + t8 * ((t56 + t203 - t2
     &15) * epinv + t319 * (ddreal(1) + epinv) + t400 * t606 - t211 + t2
     &33 + t405 - t431)) + t216 * t629 + t24 * (t8 * (t424 + t426 - t405
     & - t574 - t475 + t204 - t620) + t10 * (t384 + t341 + t223 - t472 -
     & t473)) + t3 * (epinv * (t8 * (-t56 - t203 + t215 + t400) + t10 * 
     &(-t436 + t418 - t420 + t33 + t34)) - t12 * t231 + t8 * (t211 - t23
     &3 - t400 + t431) + t10 * (t436 - t418 + t420 + t422 + t423 - t424 
     &- t426 + t297 + t268 + t245 + t30 + t31 - t33 - t34 - t204 - t626)
     &) - t46 * t9 * (t7 * (t190 * t606 - t455) + t8 * (t420 * t588 + t2
     &04 + t245 + t268 + t297 + t30 + t31 + t422 + t423 + t424 + t426 + 
     &t593 + t608 + t609 + t610 + t625)) + t10 * (-t431 - t405 + t233 + 
     &t196 - t211 - t619 + t492 - t620) + t11 - ddreal(9) * t601 * t611 
     &+ ddreal(24) * t8 * t628 + ddreal(16) * t8 * t612 - ddreal(18) * t
     &601 * t628
      t190 = t453 * t482
      t453 = t65 * t482
      t455 = t148 * t203
      t475 = t201 * t483
      t588 = (t428 + t146 + t56) * epinv
      t593 = (t233 - t199 - t209 + t571) * t12
      t601 = (-t18 * (t436 + t421) + t425 + t426) * t12
      t606 = (-t18 * (t420 + t396) + t206 + t424) * t12
      t576 = t2 * (-t622 + t616) + t24 * (t26 * (t18 * (t437 + t439) - t
     &223 - t341 - t384) + t606) - t3 * (t26 * (t18 * (-t418 - t419 - t3
     &3 - t34 - t576) + t449 * t482 + t245 + t268 + t297 + t30 + t31 + t
     &422 + t423) + t593) + t46 * (t461 + t601) - t12 * (t434 + t435 + t
     &196 - t215 + t588 + t475 + t455) - t26 * (t482 * (-epinv * t371 + 
     &t18 * (-t414 + t409 - t190) + t390 - t447 - t448 + t453) - t405 - 
     &t620)
      t29 = t114 * t461 - t129 * t26 * (-t472 + t384) + t216 * t26 * (t1
     &8 * (t33 + t437) - t223 - t30) + t24 * (t26 * (t18 * (t418 + t419)
     & + t370 + t404 - t422 - t423) + t606) - t3 * (t26 * (-t157 * t341 
     &+ t18 * t365 - t28 * t71 + t482 * (-t18 * (t454 + t414 + t408 + t4
     &10 + t193 + t195) + t375 + t377 + t378 + t390 + t440 + t449) + t23
     & - t373) + t593) + t46 * (t26 * (t18 * (-t297 + t32 + t34) - t29 -
     & t31 + t367 + t369 + t406) + t601) - t622 * t95 + t617 - t12 * (t4
     &34 + t435 + t196 + t588) - t26 * (t311 + t405) - t12 * (-t215 + t4
     &75 + t455) - t26 * (t148 * t35 - t18 * t281 + t319 * t69 + t482 * 
     &(t18 * (t407 + t409 + t379 + t371) + t482 * (-t271 + t65) - t374 -
     & t376 - t447 - t448) + (-t245 - t268) * (-t96 + t160))
      t32 = -t24 * (t384 + t341 - t472) * t26 + t3 * (t12 * (t18 * (t54 
     &+ t55) + t204 + t231 - t234 - t236) + t26 * (t18 * (t437 + t33 - t
     &477) - t223 - t297 - t30 + t618)) + t623 - t12 * ((t56 - t215) * e
     &pinv + t18 * (t159 - t400) + t196 + t203 + t211 - t431 - t432) - t
     &26 * (t590 * (t245 + t268 - t319) + t226 + t370 + t404 - t405 - t4
     &27 - t621)
      t160 = s13 + s34
      t455 = s13 - s16 + s25 + s34
      t461 = -s46 + t81
      t477 = t8 * (-t3 + t5 + t526) + t10
      t593 = s13 - s16 - s26 + s34 - s56
      t601 = s13 - s16 + s25 - s26 + s34 - s56
      t606 = t160 * t3 - t385
      t608 = s25 + s16 - s34 + s26 + s56 - t81
      t609 = -t3 * t125
      t610 = s25 + s16 - s34 + s56 + s12 + t139 - t609
      t611 = t2 * t116
      t612 = t2 * t282
      t616 = t3 * (s13 - s16 + s34)
      t617 = ddreal(9) * s25
      t618 = t79 * s13
      t619 = (s13 - t111) * s16
      t620 = t3 * t87
      t621 = -t104 * t96 - t543 + (s13 - t137 - t230) * s25 - (s12 + t12
     &1 + t117 - t81 + s26) * s26 + (t270 * t95 + t247 + t617 - t81) * s
     &34 - (-s25 - s16 - s45 - s12 - t119 - s46) * s46 + (-t3 * (s16 - s
     &46) + s45 - t115 + t304) * s56 + t184 - t618 + t619 - t620
      t531 = (t531 + t139 - t616 + s56) * s56
      t284 = (-t283 - t139 + t284 + s34) * s34
      t83 = (t83 + s16) * s16
      t622 = (t283 - t612) * s26
      t623 = (s13 - s16 - s26 + s34) * s46
      t625 = t282 * s25
      t627 = t3 * t269
      t628 = t3 * t283
      t629 = s16 - s34 + s45 + s12
      t630 = t3 * t133
      t631 = (-t176 - s16) * s16
      t632 = t3 * t104
      t187 = -t187 * t2 - (-s12 - t627 - s25) * s25 - (s16 - t135 + t628
     &) * s26 - (s13 + t258 + t111) * s34 - (t629 + t139 + t530 + s46) *
     & s46 - (-t2 * t536 - s34 + s45 + t122 + t630) * s56 - t184 - t618 
     &- t631 + t632
      t536 = s25 + s16 + s12
      t633 = t3 * t325
      t634 = s26 + s34 + s46
      t635 = t3 * t97
      t368 = (t117 + t111) * s25 + (-t24 * t270 + s13 - t121 - t252) * s
     &34 - (s25 - s13 + s16 + s45 + s12 + t119 + s46) * s46 - (t3 * t634
     & + s12 - s13 - s25 + s45 + s56) * s56 - t184 + t368 + t533 + t318 
     &+ t635
      t505 = t505 * s13
      t533 = t119 * s13
      t607 = t607 * s13
      t636 = s25 + t81
      t637 = s13 - s26 - s46
      t517 = s25 * t634 + s26 * t228 + (s16 + s26 + s56 + s12) * s34 + s
     &46 * t228 - t104 - t517
      t634 = s12 + t101
      t638 = (s25 - s16 + s26) * s25
      t639 = s12 - s34
      t640 = t162 * s46
      t641 = s12 + t84
      t270 = s16 * t270
      t642 = t3 * t58
      t643 = s25 * t186
      t644 = s13 * t67
      t103 = t10 * (-s12 * t97 + (t640 + t103) * s16 + (t24 * t87 - t3 *
     & (s25 * t149 - t140) - t334 - t93 + t537) * s26 + (-t99 * s26 - t1
     &34 * t3 + t120 + t90) * s34 + (-t132 * s56 + t128 * t2 - t3 * (t63
     &8 - t140) - s25 * t639 + t126 + t123 * s16) * s56 + t106 - t98 + t
     &108 * t634)
      t126 = s12 + s16 + s25 + s26 - s34 + s56
      t132 = s13 - s34 - s45 - s46 - s56
      t134 = s13 - s56
      t537 = t134 * t3
      t645 = epinv * t132
      t646 = s25 - s26
      t647 = s16 - s34 + s12
      t648 = (s46 + t289 + s25) * s25 + (s25 - s16 - s12 - s26) * s26 - 
     &s34 * t646 + (-t647 - t119 - s56) * s56
      t288 = (t80 + s25) * s25 + (t79 - t101) * s26 - (t249 - t119) * s3
     &4 + (t288 - t139 - s56) * s56 + t334 - t640 - t314
      t649 = s25 + s16 - s34 + s45 + s12
      t650 = t3 * t637
      t651 = (t649 - t633 + s46) * s46
      t652 = (-t176 - t101 - s25) * s25
      t653 = -(t133 + t81) * s26 + (s13 + s26 + s12 + t630 - s34) * s34 
     &+ (t180 - t650 + s56) * s56 + t184 + t618 + t631 + t651 + t652
      t654 = epinv * s34
      t655 = t165 + t402
      t656 = epinv * t416
      t83 = t468 * (t12 * ((t614 + t435) * s13 - s46 * t196 + t146 * (-e
     &pinv * t608 + s13) - t159 * (t147 * s34 + t607) - t203 * (-t590 * 
     &t67 + t607) - t215 * (-epinv * t67 + t81 * t590 + s25) - t432 * (s
     &13 + t257) + t483 * (t388 * t636 - t63 * (-t591 + t350) - t602)) +
     & t481 * (t26 * (-t23 * t653 - t268 * (epinv * t187 - t284 - t314 -
     & t531 - t618 - t622 + t623 + t625 - t83) + t311 * t648 + t405 * t2
     &88 - t443 * (epinv * t610 - t173 - t611) + t445 * t126 + t488 * t6
     &37) + t297 * t320 * (epinv * t621 + t284 + t314 + t531 + t618 + t6
     &22 - t623 - t625 + t83)) + t427 * t26 * t601 + t562 * t477)
      t147 = t249 * t451
      t284 = t586 * epinv
      t488 = t559 * t482
      t531 = t557 * t201
      t93 = t531 * (-t8 * (t216 * s25 * ((s25 + s26) * s56 + t263) - t24
     & * (-s25 * t110 + (t97 + t87) * s34 - t643 - t98 + t499 * s56 + t4
     &99 * s26) + t3 * (s12 * (t84 * t27 + s16 * (s16 - t509 - t84)) + s
     &16 * (-s25 * t67 + s46 * t67 + t104 - t186 + t270 - t362 - t498 - 
     &t644)) + t547 - t5 * ((t104 + t640 - t93) * s16 + (s25 * t641 + s2
     &6 * t133 - t334 + t620) * s26 - (t3 * t270 - t120 + t90) * s34 - (
     &(s13 + s34 - s46 - s26 - s16) * s16 - s25 * (t639 + t642)) * s56 +
     & t283 * t97 - t79 * t87 + t88 + t643 + t546)) + t103) * t482
      t83 = t468 * (t12 * (t18 * (t483 * (t557 * t517 * t508 - t132 * t2
     &2 + t147) + t217 + t587) - t233 * t593 - t431 * (-t67 + t81) - t43
     &4 * t608 + t56 * (s25 - t353) + t208 - t255 - t51 + t94 + t284 + t
     &483 * (-t219 * t637 * t481 + t218) * t69) + t26 * (-t226 * t601 + 
     &t126 * t481 * (-t495 * t245 + t488)) + t477 * t49 - t93) + t83 + t
     &468 * (t12 * (t483 * (t132 * t549 + t21 * (s25 + s46 + s26 + s45 +
     & t645 - t537) + (-t279 + t15) * (-t134 * t2 + t150 + t212) - t170 
     &* t637 * t481 - t244 * t636) + t517 * t481 * t552 * t198) + t481 *
     & (t26 * ((t281 * t648 + t288 * t319 - t35 * t653 + t415 * t637) * 
     &epinv + t126 * (epinv * t655 + t398 + t534 - t585) - t187 * t370 -
     & t441 * t610 + (-t656 - t394) * (s12 + s13 + s16 + s25 - s34 - s46
     & + s56)) + t406 * t320 * t621))
      t93 = t18 * t411
      t103 = t17 * t18
      t134 = t18 * t452
      t187 = t20 * t18
      t208 = t18 * t418
      t126 = t109 * (t26 * (t456 * (t30 * t606 + t31 * t455 + t384 * t59
     &3) + t422 - t208) + t578 * t313 * (t210 - t207) * t126 + t456 * ((
     &t236 + t424) * s13 + s34 * t235 + t160 * t433 + t199 * t253 + t430
     & * (s13 + t654) + t461 * (t206 + t426) + t483 * (t249 * (t134 - t6
     &1) - t187 * t132) + t70 * (s13 * t581 + t654)) * t12) + t468 * (t1
     &2 * (t483 * ((t93 - t450) * s13 + t132 * t16 + (t103 - t13) * (s25
     & + s46 + s26 + s45 - t537)) - t591 * (t55 + t420) - t577 * t461) +
     & t26 * (t18 * (-t33 * t606 - t34 * t455 - t439 * t593) + t481 * (t
     &28 * (epinv * t368 - (-t536 + t633 + s34) * s34 + t160 * t241 + t5
     &05 + t533) + t367 * t368 + t126 * (t538 + t397 + t401 + t497 - t27
     &5))))
      t207 = t470 * s13
      t210 = s25 + s26 + s45 + s46 + s56
      t244 = s34 - s46 - s45
      t270 = s25 + s45 + s46
      t275 = t244 * s26
      t244 = t244 * s56
      t288 = t270 * s34
      t368 = s16 * t210
      t455 = t119 * s56
      t461 = t97 + t362
      t495 = s25 * t228
      t497 = t3 * t382
      t517 = t385 - t497
      t537 = t149 + t497
      t383 = -t385 + t383
      t538 = -t191 + t84
      t549 = intHs16s25s26s34s56x1411D6eps0()
      t559 = s12 + s25 + s26
      t254 = -t254 + t84
      t578 = t3 * t191
      t581 = s25 - t578
      t585 = s12 + s13 + s16 + s26 - s45
      t593 = t24 * t382
      t601 = t3 * (s25 - s26 - s56) - t117 + t593
      t606 = -s13 + s25 + s26 + s45 + s46 + s56
      t607 = s12 + s16 + s26 - s34 - s45
      t608 = -s12 + s25 - s26 - s56
      t610 = epinv * t608
      t614 = -s25 + s16 + s26 + s56 + t111
      t621 = -s16 + s34 + s56 + t252 - t84
      t509 = t509 + t230 - t4
      t224 = t3 * t224
      t622 = -epinv * t144 - s13 + s25 + s26 + s45 + s46 + s56
      t623 = t91 + t230 - t4 + t164
      t625 = s34 + t252 - t84
      t633 = s16 - s56 - t310 + t115 - t257
      t636 = s25 + s16 - s56 - t111
      t639 = intHs16s25s26s34s56x1411D6eps1() * t18
      t640 = t466 * t456
      t643 = t456 * t463
      t648 = t26 * (t643 * t382 + t549 - t639 + t640)
      t653 = t12 * t460
      t657 = t653 * t456
      t658 = -t594 - t196
      t659 = t226 * t601
      t660 = t198 * t585
      t661 = -t18 * t22 + t1 + t14
      t190 = (-t190 + t371 + t409) * t18
      t662 = t18 * t407
      t427 = t427 * t601
      t601 = t379 * (s16 + s26 + s56 - t380 - t237)
      t105 = t487 * ((t12 * (t483 * (-t21 * (-epinv * t349 + s25 + s26 +
     & t164) + t349 * t661 + (t279 - t15) * (t212 + t151)) - t660 * t221
     & * t18) + t26 * (-t601 * t482 + t395 - t427 + t482 * (t453 - t447 
     &- t448 + t190) * t559 + (t662 - t374 - t376) * t482 * t382) + t12 
     &* (t146 + t428) * t606 * epinv + t482 * (t18 * t193 - t377) * t517
     & * t320) * propW34 * t7 + t105 * (t8 * propW34 * (t3 * (-t2 * t495
     & - t3 * t461 - t110 - t186 + t244 + t275 + t288 - t368 + t644) + t
     &5 * (t3 * (t352 + t263) + t110 + t186 - t244 - t275 - t288 + t368 
     &- t644 + t97)) - t145 * (-t110 + t97 - t186 - t368 + t644 + t288 +
     & t275 + t244 - t455)) * t167 * t483 * t482)
      t145 = t165 + t402 + t443
      t151 = t422 * t254
      t167 = t31 * t537
      t212 = t423 * t538
      t263 = t30 * t383
      t352 = t223 * t581
      t368 = (t410 + t414 + t454) * t18
      t487 = t382 * t408
      t383 = t33 * t383
      t34 = t34 * t537
      t537 = t387 * t439
      t418 = t418 * t254
      t419 = t419 * t538
      t663 = t437 * t581
      t664 = t375 * t382
      t665 = t384 * t387
      t666 = (t399 + t400 + t420) * t18
      t667 = t666 * t606
      t187 = t468 * (t12 * (t483 * (t16 * t349 + (-t103 + t13) * (t58 + 
     &t164)) - t667) + t26 * (t18 * (t482 * (t487 - t301) + t383 + t34 +
     & t537 + t418 + t419 + t663) + t482 * (-t664 + t485) - t398 - t441 
     &- t445 - t665 - (t390 + t440 + t449) * t482 * t559)) + t468 * (t12
     & * (-t144 * t231 + t199 * t622 + t606 * (t209 + t211 + t424) - t18
     &7 * t349 * t483) + t26 * (-t28 * (-epinv * t621 + t497) - t341 * (
     &-epinv * t625 - t578) - t365 * (s12 + t610 + t84) + t367 * t621 + 
     &t369 * t625 - t373 * t608 - t69 * t145 - t151 - t167 - t212 - t263
     & - t352) + t368 * t26 * t482 * t559)
      t105 = t166 * t470 * t606 * (t577 - t206 - t426) + t169 * t109 * (
     &t657 * t606 + t648) - t174 * t187 + t197 * (t105 + t468 * (t12 * (
     &t144 * (t233 + t431) + t606 * (t434 + t435) + t607 * t658 + t622 *
     & (-t203 - t215) + t660 * (t18 * t222 + t175 - t189)) + t26 * (-t23
     & * t636 + t245 * (-epinv * t509 + t224 + t252) + t268 * (-epinv * 
     &t623 + t229 + t252) - t281 * (s25 + s12 + t610) + t297 * (epinv * 
     &t633 - t386 + t593) - t311 * t608 + t319 * (epinv * t614 - s12 + s
     &16 + s26 - s34 + s56) - t35 * (epinv * t636 + s12 + s34) - t370 * 
     &t623 - t404 * t509 + t405 * t614 + t406 * t633 + t659)))
      t109 = s12 + s34 - s56
      t144 = s12 + s13 + s16 + s26
      t187 = t176 * t3 + s16 - s25 + s26 + s34 - s56
      t301 = -s12 - s13 + s25 - s34 + s56
      t349 = t101 * s12 + (t80 + t81 - s25) * s25 + (t2 * t80 - s25 + s2
     &6 + t81) * s26 + (-s16 + s26 - s12 + t530 - s34) * s34 + s46 * t58
     &5 + (-t3 * (s25 - s13 - s34) + s16 - s26 - s45 + s12 - s56) * s56 
     &+ t182 - t183 - t184 + t87
      t485 = t3 * t79
      t509 = (-s16 - t121 + t81) * t149
      t321 = s12 * (s12 - t321)
      t530 = t642 * t149
      t585 = t149 * t58
      t608 = s12 * t162
      t610 = t2 * t58
      t614 = (-t250 - t81) * t97
      t621 = s12 * t104
      t622 = s12 * t87
      t623 = t2 * t646
      t625 = t216 * t97
      t136 = -propZ25 * (t3 * (t513 * propW16 * s16 * t67 + t98) + (s25 
     &* (-s16 + t81) - t265) * s16 + ((t248 + t111) * s26 + t111 * t133 
     &- t265 - t509) * s26 - (t321 + t530) * s34 + (t608 + t585) * s46 +
     & ((t253 - t99 - t139 + t257 - s56) * s56 - t3 * (-(-t79 - t101) * 
     &s25 + t110 + t87) + t136 + t543 + (s13 - t545 + t115 + t111) * s26
     & + (-s25 - s13 - s12 + t99 + t139 - s34) * s34 + (s25 + s34 - s26 
     &+ s12 - t101) * s46 - t265) * s56 + t614 + t621 + t622) * t2 + t62
     &5
      t143 = t10 * (-(s25 * (s16 - t485) - t265) * s16 - (-(t248 - t111)
     & * s26 + (-s13 + t99 + s12) * s12 + t509) * s26 + (t321 - t530) * 
     &s34 + (-t608 + t585) * s46 - ((-t125 - t610 + t257 - s56) * s56 - 
     &t3 * (-s25 * t269 + t110) - t543 + (t3 * t80 + s13 - t86) * s26 + 
     &(-t176 + t610 - s34) * s34 - (-s12 + s25 + s26 - s34) * s46 - t265
     & - t619) * s56 + t614 - t621 - t622 + t143)
      t513 = s12 + s16 + s26 + s46
      t545 = epinv * t270
      t500 = t3 * t501 + t139 + t500
      t501 = s12 + t81
      t619 = -t104 - t265
      t633 = t2 * t176
      t636 = (s13 + t230) * s16
      t185 = -t185 * t2 - t24 * t265 + t314 + (t633 - t101) * s25 + (-s1
     &3 - t4 + t278) * s26 + (-t2 * (s13 - s26 + s12) + t101 + t86) * s3
     &4 + (-s25 + s34 + s26 + t230) * s46 + (t125 + t247 + t139 - t429 +
     & s56) * s56 + t636
      t278 = t46 * s13
      t660 = t540 * s25
      t131 = t46 * t619 - t315 + (t247 + t135 - s16) * s16 + (t307 - t13
     &7 + t252) * s25 + (s12 * t540 + s26 - t117 - t131) * s26 + (t660 +
     & t259 - t278 + t113 - t230) * s34 + (-t272 + t310 - t4 + t257) * s
     &46 + (s12 * t216 + s46 * t3 + s56 + t139 - t269 + t304 - t86) * s5
     &6
      t272 = t3 * t186
      t287 = (s16 + t252 - t81) * s26 - s34 * t287 - (s25 - s26 - s12) *
     & s46 + (-t3 * (s13 + s34 - s12) + s16 + s25 + s46 + t123) * s56 + 
     &t290 + t291 + t272 + t532
      t246 = t3 * (t186 - t265) + t292 - (-t246 + s25) * s25 + (s25 + t3
     &10 + t117 - t81) * s26 - (-t149 + t330) * s34 + (t127 + t442) * s4
     &6 + (-t3 * t344 + s25 + t173 + t230) * s56 + t87 + t252 * s16
      t291 = t2 * s46
      t129 = t110 * t96 + t24 * (-t97 + t186 - t265) - t318 + (t112 - t1
     &17 + t111) * s25 + (s12 * t129 + s16 - t391) * s26 + (s26 + t100 -
     & t633 + t101) * s34 + (-t3 * t646 + s34 + t230) * s46 + (-s34 + t1
     &42 + t312 - t112 + t291 + t101) * s56 + t636
      t142 = t2 * (s13 - s16 - s46 - s56) - t173 - t429
      t304 = t178 * s46
      t179 = s13 * (t385 + t179)
      t312 = t3 * (s16 * t178 + t362)
      t315 = t2 * (-t97 + t182 - t104)
      t330 = t24 * t361
      t344 = t96 * s26 * t178
      t361 = t499 * t46
      t264 = -t228 * (-t264 - t139 - t502) + t366
      t362 = t24 * t127
      t366 = t501 * s16
      t391 = s16 + s26 - s34 - s45 - s46 + s56
      t7 = t465 * propW34 * t7
      t465 = t574 - t426
      t502 = t571 - t209
      t540 = t569 - t206
      t103 = -t199 * (s25 - s56 + s12 - t545) + t231 * t270 + t476 * t10
     &9 + t559 * t502 + t540 * (t522 + t628) + t465 * (-t3 * t446 + t227
     &) + t103 - t16 * s34 * t483
      t285 = t551 * (t341 * (epinv * (-s12 * t142 - t110 + t179 - t186 -
     & t260 - t261 - t304 - t312 + t315 - t330 - t344 + t361) - t3 * (t1
     &82 - t87 + t340 + t339 - t336 + t337 - t338 - t335 + t334) + t292)
     & + t365 * (epinv * t264 + s25 * t285 - t293 * t641 - t3 * (s26 * t
     &641 + t97) - t290))
      t91 = t468 * (t12 * t103 + t26 * (t481 * (t28 * (epinv * ((-s34 * 
     &t95 + t173) * s25 + t2 * (s25 * s56 + t104) + t24 * (t97 - t498) +
     & t3 * (-s56 * t91 + t358) + s13 * (s56 - t115 + t302 - t101) - (-t
     &2 * t516 - s16 - t115 - t241) * s16 - t359 - t360) + (t3 * (s25 + 
     &s46 + t239 + t654) + t354 + t355 - t356 - s34 + t351 - t350 + t363
     &) * s12 - t364) + t367 * (-t2 * t619 + t543 - (t24 * t282 + s12) *
     & s25 + (t286 * t46 + t117) * s26 + (-t100 - t168 + t612 + t111) * 
     &s34 + (t127 - t252 + t84) * s46 - (t3 * t91 - s13 - s16 + s46 + s5
     &6 + t230 - t4) * s56 + t87 - t101 * t176) + t369 * (s12 * t142 + t
     &110 - t179 + t186 + t260 + t261 + t304 + t312 - t315 + t330 + t344
     & - t361) - t373 * t264) + t482 * (t144 * t440 - t301 * t449) + (-t
     &595 + t31) * (s25 - s16 + s56 + t497)) - t285)
      t103 = t18 * t408
      t142 = t18 * t195
      t129 = -t245 * (-t3 * ((t249 - t117 + t111) * s26 - t104 - t186 - 
     &t87 + t97) - t267 + t292 + epinv * t185 + (-t247 + t81) * s25 + (s
     &12 - t173 - t99 + t81) * s34 + (-t3 * (s25 - s16 + s34 - s26) - t2
     &52) * s46 - (s12 - t173 - t362 + t609) * s56 - t366) - t268 * (t3 
     &* (t110 + t104 + t87) - t267 + epinv * t129 - (-t269 + t230 + s25)
     & * s25 - (-t250 * t96 + t249) * s26 + (-t286 - t123 - t99 + t81) *
     & s34 + (-t646 - t252 + t322) * s46 - (-t162 * t2 - s56 + t111 - t1
     &25) * s56 - t366) - t297 * (epinv * t131 - s25 * t298 - s34 * t300
     & - s46 * t303 - s56 * t305 + t316 - t317 + t318) + t311 * t287 + t
     &319 * (epinv * t246 - t182 + t314 - t334 + t335 + t336 - t337 + t3
     &38 - t339 - t340 + t87) - t370 * t129 - t404 * t185 + t405 * t246 
     &- t406 * t131 + t655 * (epinv * t513 - t173 - t611)
      t131 = t146 * epinv
      t179 = t492 - t211
      t65 = t468 * (t12 * (t109 * (epinv * t428 + t435) - t203 * (epinv 
     &* t391 + s12 + s25 + s26) - t215 * (s25 + s26 + s12 - t545) - t233
     & * t391 + t431 * t270 + t434 * t559 + t483 * (s34 * t661 + t279 * 
     &t515) + t179 * (t523 + t111) + t198 * t66 * t481 * t349) + t26 * (
     &t481 * (-t23 * ((-t230 + t258) * s26 - t3 * (-t97 - t265) + (-t501
     & + s16) * s16 + (t325 - t4 - t251 + s34) * s34 + (t133 - t111) * s
     &46 - (s34 - t326) * s56 - t273) + t395 * t503) + t482 * (t371 * (e
     &pinv * t301 + s12 + s13 + s16 + s26) + t144 * t482 * (-t271 + t65)
     &) + t659) + t390 * t320 * t187 * t482)
      t109 = -t374 - t376
      t136 = t468 * (t12 * (t481 * (t553 * t333 * t483 + t323 * t349) + 
     &t53 * t658 - t555 - t556) + t26 * (t482 * (s12 * t109 + t69 * (-t2
     &20 * t503 * t481 - t68) - t490 - t491 - t78 - t561 * t301) - t427 
     &- t392 * t333 * t481 + t550 * t481) + t527 + t563 - t531 * (t8 * (
     &t136 * s25 + t24 * (t110 * (s12 + t248) + t614) + t3 * ((-s25 * (s
     &12 - t112 + t101) - t265) * s16 + (-t3 * t509 + s12 * (s12 - t282 
     &+ t4)) * s26 - (t585 * t24 + t321) * s34 + (t608 + t530) * s46 + (
     &(-t125 - t137 + t623 + t257 - s56) * s56 - (s12 - t539) * s16 + (-
     &t130 - t112 + t252) * s25 + (s13 - t541 + t617 + t111) * s26 + (-t
     &176 + t137 - t623 - s34) * s34 + (s34 - s26 + s12 - t117 + t84) * 
     &s46 - t265 - t544 + t625 - t314) * s56 + t621 + t622) + t547) + t1
     &43) * t482)
      t65 = t468 * (t12 * (t558 * t170 + t192 - t47) + t26 * (t481 * (ep
     &inv * (t118 * t514 - t333 * t415) - t535) + t482 * (t194 + t274) *
     & t18) - t564 * t320 * t187 * t482) + t136 + t65 + t468 * (t12 * (-
     &t483 * (t15 * t515 + t21 * (s25 + s26 + s56 - t654)) + t131 * t559
     &) + t481 * (t129 * t320 + t26 * (-t281 * (-epinv * t287 + t290 + t
     &294 + t295 + t296) - t35 * (epinv * (t3 * s25 * t200 + t2 * s25 * 
     &t516 + s13 * (s34 - t630) + (-t3 * t494 + s16 + s46 + t610) * s16 
     &- s34 * t348 + t358) + (t3 * (-epinv * t293 + s26) - s34 + s46 + s
     &56 + t238 + t239 + t240 - t242 + t243) * s12) + t441 * t500 + t443
     & * (epinv * t500 - t259 - t444) + t513 * (t398 + t445))) + t26 * (
     &t447 + t448) * t482 * t301)
      t13 = -t169 * t7 * (t464 * (t26 * (t643 * t381 + t549 - t639 + t64
     &0) + t657 * (s12 + s25 - s56)) + t456 * ((-t565 - t49 + t524 - t56
     &2) * s16 + t529 + t566) * t519 * propW16) - t174 * (t468 * (t12 * 
     &(t237 * t20 * t483 - t13 - t50) + t26 * (t254 * (-t208 + t422) + t
     &474 * t581 + t482 * ((-t103 + t375) * s12 + t18 * (-t144 * t454 + 
     &t301 * t410) + t191 * (t142 - t378)) + t538 * (-t626 + t423) + t47
     &9 + t480) + t26 * (t656 + t394) * t481 * t333) + t91) - t197 * t65
     & - t518
      t65 = t116 + t139
      t91 = s13 - s45 - s46 - s56
      t118 = s13 + s16 + s26 - s34 - s45 - s46
      t129 = t3 * (s13 - s34 - s46 - s45) + s16 + s26 - s56
      t136 = s16 + s26 - s34 + s46 + s56
      t116 = t116 + t168
      t137 = propW16 * s16
      t143 = t137 * t519
      t144 = t8 * (t2 * (t143 + t9) - t3) + t10
      t182 = t515 + t135
      t185 = t3 * t600
      t187 = -t389 + t539 + t185
      t194 = s25 - s13 + s16 - s34 + s46 + s56 + s12 + t507 + t119
      t200 = t2 * t75
      t238 = t539 - t200
      t77 = t135 - t77
      t239 = t3 * t649
      t240 = s46 + s26
      t241 = t114 * s13
      t242 = (s46)**2
      t243 = (s13 + t117 + t111) * s25
      t246 = (t501 + s16) * s16
      t130 = -t186 * t96 + t3 * (t104 + t246 + t243) - t46 * (t242 + t18
     &4 + t618) + t532 + t543 + (s13 * t306 + s12 + t130 + t660) * s26 -
     & (t24 * t269 + t111 + t181 + t86) * s34 - (s26 * t306 + t180 * t46
     & - t241 + t362 + t84) * s46 - (s45 * t46 + t125 * t306 + t127 * t2
     & + s25 + t247 + t259) * s56
      t177 = t177 * s13
      t247 = t242 + t186 + t184
      t250 = t24 * t247
      t251 = t2 * t629
      t260 = t95 * t325
      t261 = t2 * s45
      t264 = t95 * s46
      t265 = t24 * t180
      t267 = (t99 + t252) * s25
      t259 = -t2 * (-t97 + t618) - t250 + t632 - (-t342 - t278 - t117 + 
     &s26) * s26 - (s13 + t173 + t121 + t117 + t111) * s34 - (s25 + t259
     & - t307 + t265 + t117 - t257) * s46 - (s46 * t216 + s25 + t117 + t
     &181 - t257 + t265 - t307) * s56 + t246 + t267
      t269 = t289 * s13
      t270 = (t80 - t539 + s25) * s25
      t273 = t242 + t184
      t274 = -s25 + t112
      t122 = -t2 * (-t97 + t186 + t618) - t24 * t273 + t3 * t343 + (s12 
     &+ t617 + t539 + t99) * s26 - (s13 + t122 + t121 + t117 + t111) * s
     &34 - (s25 - t260 + t265 + t117 - t257) * s46 + (-s45 * t24 + s34 -
     & t101 - t123 - t252 - t264 + t278) * s56 + t246 + t267
      t265 = t2 * t180
      t267 = -t2 * (-(s16 + t589) * s26 + t184 + t186 + t242 + t618) + t
     &635 + (-s12 - t168 - t4 - t627 + s34) * s34 - (-t325 * t46 + s25 +
     & t265 + t322) * s46 - (-t46 * t637 + s25 + t265 + t322) * s56 + t2
     &43 + t246
      t282 = (t649 - t650 + s56) * s56
      t285 = t112 * s26
      t286 = s25 - s13 + s26
      t272 = -(t536 - t135) * s26 - s34 * t325 - (s46 + t649 + t139 - t8
     &1) * s46 - (t2 * t125 + t3 * (s25 + s16 - s34 + s12) + s45 + t173)
     & * s56 - t184 + t505 - t272 - t532
      t82 = epinv * t272 + (t82 + s25) * s25 - s34 * t299 + (t286 * t3 +
     & s46 + t629) * s46 + (t3 * (s25 - s13 + s46 + s26) + t629 + s56) *
     & s56 + t184 - t269 + t119 * t299
      t180 = t3 * t180
      t258 = (t112 + t258) * s26 - t3 * (t242 + t186 + t184 + t618) + (-
     &t176 - t139 - t630 + s34) * s34 - (t141 - t506 + t180) * s46 - (-t
     &24 * t637 + t141 + t180) * s56 - t631 - t652
      t99 = -t2 * (-t110 + t242 + t184 + t618) + t3 * (t97 - t186) + (s1
     &2 + t100 + t539 + t99) * s26 + (-s12 - t123 - t4 - t627 + s34) * s
     &34 - (s25 - t278 + t123 + t265 + t322) * s46 - (t96 * t125 + t127 
     &+ t261 + t442) * s56 + t243 + t246
      t100 = (t240 * t3 + s56 - t135 + t649) * s56
      t123 = (t649 - t135 + t119 + s46) * s46
      t140 = t3 * (t110 - t242 - t184) - (-t634 - s25) * s25 + (t133 * t
     &24 + t176) * s26 - (s12 + t168 + t630 - s34) * s34 - (-t2 * t325 +
     & t141 + t180) * s46 - (-t3 * (s13 - s45) + t291 + t329 + s56) * s5
     &6 + t140 - t183 + t335 + t87
      t176 = t79 * t133
      t139 = (t139 + t101) * s16
      t180 = s16 * (t516 + s16)
      t183 = (-t289 - s25) * t97
      t98 = t10 * (s25 * t87 + ((s16 + t4) * s26 + t108 - t176 + t543 + 
     &t544) * s26 - ((s25 + s16 + s26) * s25 - t120 + t139) * s34 + (t90
     & + t180) * s46 + (s56 * t133 + t3 * ((s25 + s16 - s34 + s26) * s16
     & + t97) + s46 * t133 - t176 - t499 + t115 * s26) * s56 + t97 * t28
     &9 + t106 + t107 + t98)
      t107 = t136 + t112
      t108 = -epinv * t95 + t24
      t133 = -t24 * t75 + t539
      t176 = t112 - t200
      t58 = t26 * ((t441 + t445) * s26 - t130 * t406 + t194 * t402 + t23
     & * t258 - t245 * (t2 * t247 - t24 * t269 - t3 * t274 * s26 + epinv
     & * t122 + s34 * t274 + (-t307 + t173 + t115 + t251) * s46 + (t240 
     &* t46 + t115 + t251 - t307) * s56 + t270) - t268 * (t2 * t273 + t2
     &4 * (t186 - t269) + t314 + epinv * t259 + (-ddreal(9) * s13 + t4 +
     & t80) * s26 + (-t58 + t112) * s34 + (-t260 + t115 + t251) * s46 + 
     &(t24 * t647 + t121 + t181 + t261 + t264 - t309) * s56 + t270) + t2
     &81 * t82 - t297 * (epinv * t130 - t46 * t177 + t3 * ((-t539 + t168
     & + t239) * s46 + (t24 * t240 + t239 - t539) * s56) + t250 - t241 *
     & s26) + t311 * t272 + t319 * (epinv * t140 - t3 * t177 + t100 + t1
     &23 + t184 - t285) + t35 * (epinv * t258 - t177 + t184 + t282 - t53
     &3 + t651) - t370 * t259 - t404 * t122 + t405 * t140)
      t58 = t468 * (t12 * (-t107 * t233 - t146 * (-epinv * s16 - t81 * e
     &pinv + s13 + s25 + s26 + s56) + t159 * ((t24 - t205) * s13 + s25 +
     & s26 + s56 + t654 + t257) + t187 * t211 + t196 * (s46 + t135) + t2
     &03 * (-epinv * t107 + s25 - s26 - s56 + t112 - t322) + t215 * (s13
     & * t108 + t115 + t345 - t353) - t431 * (s46 + t307) - t432 * (-s34
     & + t539) + t434 * (s16 + t81) - t56 * (-t135 * epinv + s13 + s25 -
     & t353)) + t26 * (t133 * t226 + t481 * (t116 * (-t394 - t401) + t19
     &4 * t443 + t65 * (t392 - t395 + t398) + t534 * s26) + t482 * (-t12
     &9 * t193 + t136 * t414) * t18 + t377 * t129 * t482)) + t468 * (-t1
     &8 * (t400 * t12 * t187 + t438 * t26 * t133) + t26 * (t109 * t132 +
     & t136 * (-t371 - t390) + t379 * (s13 + s16 - s34 - s46 + s26 - s45
     & - t645)) * t482 + t58 * t481)
      t79 = t531 * (t8 * (-t24 * (t110 * t102 - t183) + t3 * (((s25 + s1
     &6 - t79) * s16 + t104) * s16 + ((s16 + t485) * s25 - t120 + t544 -
     & t625) * s26 - (-(-s25 + s13 - s12) * s16 - t138 + t117 * s26) * s
     &34 + (t180 - t92) * s46 + (s56 * t85 + (-t342 + t113) * s26 + t3 *
     & (t87 - t512 + t124) - t543 + s12 * t85 - s16 * t253) * s56) + t54
     &7 + t5 * ((-t120 - t638 + t139) * s34 - (-t90 + t180) * s46 - (-t2
     &4 * s26 * t149 - t3 * (t97 - t128) - s46 * t149 + t149 * t79 + t49
     &9) * s56 - t149 * (-s26 * (-t79 + t115 + t117) - t186) - t87 * (s2
     &5 - t79) - t106 - t88 - t183 - t546 + t248 * t110)) - t98) * t482
      t80 = t70 * (t27 * t3 + t347)
      t85 = t482 * t18
      t87 = t468 * (t12 * (-t182 * t234 - t199 * (t596 - t353 + t84) + t
     &393 * t540 + t483 * (t249 * (-t134 + t61) + t346 * (-t93 + t450)) 
     &+ t502 * t580 + t496) + t144 * (t562 + t48 + t49) + t154 * (-t143 
     &* t2 * t8 + epinv * t144) + t26 * (-t176 * t30 - t223 * t77 - t238
     & * t384 + t91 * (-t422 - t423) + t85 * (t118 * t195 + t132 * t408 
     &- t136 * t454) + t313 * t417 * t65))
      t33 = t87 + t468 * (t12 * (t160 * t476 + t231 * (s46 + t112) - t43
     &0 * t592 + t465 * (t227 + t81) + (t478 - t236) * (s34 + t81) - t80
     & + t469 * t182) + t26 * (t18 * (t176 * t33 + t238 * t439 + t437 * 
     &t77 + t484 * t91) + t481 * (t267 * t367 + t28 * (epinv * t267 + t3
     & * (-t177 + t282 + t184 + t651) - t285) + t341 * (epinv * t99 - t2
     &4 * t177 + t3 * (t123 + t100 + t184) - t309 * s26) + t365 * t82 + 
     &t369 * t99 + t373 * t272 + t397 * t65) + t482 * (-t118 * t378 - t1
     &32 * t375 + t136 * t440)))
      t77 = s13 * t433
      t82 = ((t565 - t524) * s16 - t529 - t566) * t519 * propW16
      t33 = t166 * t468 * t26 * t91 * (-t31 + t595) + t169 * t568 * (t46
     &4 * (t12 * (t521 + t77 + t510) + t26 * (t36 * t91 + t458 * (-t75 +
     & t81))) + t82) + t174 * t33 - t188 * t470 * (t50 + t51) - t197 * (
     &t468 * (t12 * (-t483 * (t280 + t158 + t542) - t192 + t255 + t47 - 
     &t94) + t662 * t26 * t132 * t482 + t412 * t12 * t65 * t332) + t468 
     &* (t12 * ((t153 - t586) * epinv + t157 * (t156 + t217) - t18 * t58
     &7 + t582 - t605 - t603 * t483 + t604 * t483 + t602 * t483 + t557 *
     & (t554 * t66 - t170 - t553) * t65) + t481 * (epinv * (t116 * (-t17
     &2 - t416) + t65 * (t165 + t415)) + t488 * s26) * t26 - t79) + t58)
     & - t518
      t36 = epinv * t70
      t1 = -t166 * t468 * (t12 * (t483 * (t18 * (-t20 + t411 + t413 + t4
     &52) + t16 - t388 - t450 - t61) + t199 + t203 + t209 + t211 + t424 
     &+ t430 + t433 - t666 + t36) + t26 * (-t262 + t28 + t30 + t31)) + t
     &169 * t470 * t483 * (t279 - t15) + t174 * t468 * (t12 * (t483 * (t
     &18 * (-t22 + t63 + t451) + t21 * t69 + t1 + t14 - t59 - t62) + t59
     &0 * (-t159 - t215 + t428) - t431 - t432 + t434 + t435 + t131) + (t
     &44 + t23) * t26) + t372 * t468 * (t653 + t567) + t89 * t470 * (t18
     & * (t17 * t483 + t396 + t436) - t206 - t426 - t613)
      t14 = (-s25 - s34 - s45 - s46 + s16) * s16
      t15 = t72 * t162
      t16 = s46 + s26 + s56 + s45
      t17 = t2 * t228
      t20 = s16 - t252 + t429
      t21 = t178 - t111
      t22 = t141 - t111
      t30 = t607 + t81
      t31 = t583 * s12
      t44 = t389 - t230 + t429
      t58 = s12 + s26 + s56
      t65 = t228 - t310 + t357 + t101
      t66 = t385 - t252 + t429
      t75 = s25 + s16 - s34 + s46 + s56 + s12
      t79 = t3 * (s13 + s26) + t75
      t87 = t75 + t135 + t119
      t75 = t75 + t539 + t119
      t88 = t210 - t81
      t90 = s25 - s16 + s46 + s26 + s56 + s45 - t135
      t91 = epinv * t58
      t92 = -t217 - t587
      t93 = t483 * t18
      t98 = t377 * t320 * t517
      t14 = t475 * (t8 * (-t24 * t461 - t3 * (t110 + t186 - t14 - t644 -
     & t288 - t275 - t244) - t46 * (-t137 * t27 * t67 + t495) + t5 * ((s
     &25 - s13 - s16) * s16 + (-s13 + t84 + s26) * s26 - (s25 - s13 - s1
     &6 + s26) * s34 + (-s13 - s34 + s46 + s45 + t642 + s56) * s56 + t97
     & + t15 - t101 * propW16 * t27 * t67)) - t10 * (-t110 + t97 - t186 
     &+ t14 + t644 + t288 + t275 + t244 - t455)) - t98
      t14 = t468 * (t12 * (-t196 * t30 + t520 * t92 + t56 * (-epinv * t3
     &0 + s13 + s25) + t93 * (-t413 * t599 + t554 * t508)) + t26 * (-t23
     & * t22 + t245 * (t3 * (s25 * t590 - s26 - s56 + t237) + s16 * t69 
     &+ t31) + t319 * (-epinv * t21 - s12 + s16 + s26 - s34 + s56) - t35
     & * (epinv * t22 + s12 + s34) + t370 * t44 + t404 * (s16 - t230 + t
     &429) - t405 * t21 + t482 * (t559 * (t453 + t190) - t601 + t662 * t
     &382) + t659) + t482 * t14)
      t21 = t109 * t382
      t22 = t468 * (t12 * (t483 * (t249 * (-t59 - t62) + t388 * t599 + t
     &69 * (t218 + t219) + t158 + t170 + t603) + t284 - t556) + t26 * ((
     &t172 - t415 + t416 + t417) * epinv + t482 * (t489 * t559 + t21) - 
     &t392 + t394 - t427) + t477 * (t562 + t49))
      t30 = t193 * t320 * t517
      t14 = t468 * (t12 * t94 + t18 * (t147 * t12 * t483 + t30 * t482) +
     & t26 * (t395 + t397 + t398 + t401 + t441 + t445)) + t22 + t14 + t4
     &68 * (t12 * (t146 * (epinv * t90 + s13 + s25 + s26 + s56) - t159 *
     & (t3 * t515 + t584) + t179 * (t228 * t24 + t3 * (s13 + s34 + s46 +
     & s45) - s16 + t86) + t203 * (epinv * t87 - s26 - s45 - s46 - s56 -
     & t4 - t616) + t215 * (-t2 * (s13 + s26 + s56) - t121 + epinv * t75
     & - s45 - s46) + t233 * t87 + t428 * (epinv * t88 + s13 + s34) + t4
     &31 * t75 + t434 * t90 + t435 * t88) + t154 * (epinv * t477 - t528)
     & + t26 * (t145 * t403 + t268 * (epinv * t44 + t229 + t252) - t281 
     &* (s25 + s12 - t91) + t297 * (epinv * t65 - t386 + t593) + t311 * 
     &t58 + t406 * t65) + t477 * t48 + t12 * t552 * t198 * t332)
      t22 = t468 * (t12 * (t93 * t276 + t255 + t459 + t50 + t51) + (t472
     & - t384) * t387 * t26 + t26 * (-t142 + t378) * t482 * t191)
      t4 = t22 + t468 * (t12 * (t511 * (-t478 + t236) + t25 - t469 * t57
     &9) + t26 * (t18 * (t34 + t383) + t254 * (t208 - t422) + t482 * (t3
     &82 * (t103 - t375) + t559 * (t368 - t390 - t440 - t449)) + t538 * 
     &(t626 - t423) + t581 * (t473 - t223))) + t468 * (t12 * (-s13 * t23
     &5 - t199 * (epinv * t79 - s13 - s26 - s45 - s46 - s56 - t4) - t231
     & * t79 - t430 * (-t591 - t257) + t570 * (s13 + t121 + t17 + t150) 
     &+ t572 * (t72 + t115 + t17 + t81) + t573 * (-s16 + s46 + s26 + s56
     & + s45 + t185) + t575 * (t16 * t3 + t115 - t127) + t80) + t26 * (-
     &t20 * t367 - t28 * (epinv * t20 + t497) - t341 * (epinv * t66 - t5
     &78) - t365 * (s12 - t91 + t84) - t369 * t66 + t373 * t58 - t167 - 
     &t263))
      t4 = t169 * t7 * (t464 * (t456 * (t460 * (t16 + t84) + t510 + t560
     & + t77) * t12 + t648) + t82 * t456) - t174 * t4 + t197 * t14 - t51
     &8 + t207 * t166 * t432
      t14 = (t396 + t436) * t18 - t206 - t426
      t16 = t3 * ((t399 + t400 + t420) * t18 - t199 - t209 - t211 - t424
     &)
      t17 = t24 * t14
      t20 = -t46 * t327
      t22 = t2 - t467
      t25 = t24 * t456
      t44 = t22 * propZ25 - t25
      t25 = -t22 * propZ25 + t25
      t58 = t2 + t467
      t65 = t3 * t456
      t66 = t58 * propZ25 - t65
      t75 = -s13 + s25 + s26 - s34 + s45 + s46 + s56
      t58 = -t58 * propZ25 + t65
      t65 = propZ25 * (-t2 - t467 - t525) + t456 * (t3 + t525)
      t77 = s13 - s45 - s46
      t79 = s12 + s16 + s25 + s26 - s34 - s45 - s46
      t80 = propW16 * s34 + ddreal(1)
      t81 = -t348 * t210
      t82 = -s13 + s16 - s34 + s56 + s12 + t642
      t86 = -s13 + s16 - s34 + s12 + t610 + t164
      t87 = -s13 + s56 + s12 + t642
      t78 = t486 * (t12 * (-t156 * t18 - t196 * t79 - t203 * (-epinv * t
     &87 - s13 + s25 + s26 + s45 + s46 + s56) - t215 * (-epinv * t86 - s
     &13 + s25 + s26 + s45 + s46 + s56) + t233 * t87 + t431 * t86 + t432
     & * t515 - t434 * t77 + t435 * t75) + t26 * (t245 * (t31 + t224) - 
     &t281 * (s25 + t380) + t405 * (s16 + t111) + t482 * (t78 + t171) + 
     &t520 * (t165 + t402 + t443) + t659))
      t21 = t486 * (t12 * ((-t146 * t77 + t159 * t515 + t428 * t75 - t56
     & * t79) * epinv + t18 * t92) + t26 * (t268 * (t31 + t229) + t319 *
     & (s16 * t590 + s26 - s34 + s56 + t152) - t35 * (s34 - t152) + t482
     & * (t225 * t390 + t283 * t489 - t371 * (s16 - s34 + s26 + s56 - t3
     &80 + t350) + t69 * (t68 + t220) + t21 - t601 - t453 * t191) - t427
     &) - t98 * t482)
      t5 = t21 + t58 * t94 + t65 * (t48 + t49) + epinv * (t155 * t65 + t
     &58 * t586) + t154 * (epinv * t65 - t526 * t456) + t44 * ((t172 + t
     &416 + t417) * epinv + t394 + t397 + t401) + t66 * (t255 + t50 + t5
     &1 + t47 + t555) + t311 * s12 * t25 + t486 * (t18 * (t26 * (-t225 *
     & t414 + t283 * t409 + t382 * t407) + t30) + t475 * (t8 * (t3 * (-s
     &16 * (t2 * t27 + t72) + t644 + t81) - t5 * (s16 * (-t3 * (-propW16
     & * t110 - propW16 * t186 + s25 * (-propW16 * t348 + ddreal(1)) + s
     &26 * (-t164 * propW16 + t80) + s56 * t80) - t72) + t644 + t81 + t6
     &20 * propW16 * t27) + t525 * t27 * t67) + t10 * ((t299 + s26) * s2
     &6 - s34 * t286 + (s25 - s13 - s34 + s46 + s45 + t119 + s56) * s56 
     &+ t15 - t334))) * t482 + t271 * t486 * t26 * t191 * (t482)**2 + t7
     &8
      t8 = t44 * s12
      t10 = t486 * (t12 * (-t199 * (-epinv * t82 - s13 + s25 + s26 + s45
     & + s46 + s56) + t231 * t82 + t667 + t36 * t27) + t26 * (t18 * (-t6
     &63 - t383 - t34 - t537 - t418 - t419) - t28 * (t76 - t257) + t365 
     &* (t380 + t84) + t482 * (-t191 * (t378 + t440) + t283 * t449 + t66
     &4) + t151 + t167)) + t25 * (s12 * t23 + t398 + t441 + t445) + (epi
     &nv * t430 + t235) * t66 * s34 + t8 * t373 + t486 * (t12 * (t27 * t
     &433 - t606 * (t209 + t211 + t424)) + t26 * (-t297 * (s16 + s26 + s
     &56 + t76 - t257) - t341 * (t3 * t67 + t76) + t212 + t263 + t352 + 
     &t665) + t85 * t26 * (t191 * (t454 + t195) - t283 * t410 - t487))
      t15 = struc5PP + struc6PP
      t2 = (t483 * (t483 * (t18 * (t221 - t222) - t175 + t189) - t201) -
     & t196 + t203 + t215 - t434 - t435 - t588 - t20 - t17 - t16) * stru
     &c4PP - (-t3 * (t231 + t204) - (t203 + t215 - t56) * epinv + t196 -
     & t233 - t431) * struc7PP + (-t2 * (t196 + t594) + t216 * t231 + t2
     &33 * t96 + t24 * (t148 * t199 + t18 * (t436 + t400 + t54) - t211 -
     & t234 - t426 + t432) - t3 * (t18 * (t430 - t420 - t55) + t483 * (t
     &18 * (t452 + t413 + t411) - t388 - t450 - t61) - t70 * t71 + t131 
     &- t235 + t236 + t424 + t434) + t431 * t95 + t46 * (t18 * (t421 + t
     &399 + t396 + t52) - t206 - t209 - t232 - t425 + t433) + t159 * t16
     &1 + t18 * t428 + t203 * (-t2 + t205) - t215 * t108 + t483 * (t18 *
     & t63 + t277 - t59 - t62) + t146 - t435 + t56) * struc9PP + t15 * (
     &-t434 - t435 - t196 + t203 + t215 - t588 - t475 - t20 - t17 - t16)
      t2 = t468 * (t12 * t2 + (t3 * (t12 * (t231 + t213 + t204) + t214) 
     &- t12 * ((t56 - t146 - t428) * epinv + t196 - t434 - t435 + t42) -
     & t202 - t266 - t308) * struc11PP + (t3 * (t12 * (t18 * (t436 + t39
     &6 - t199) - t206 + t213 + t231 - t426) + t214) - t266 - t12 * (t42
     & - t434 - t435 + t324) - t202 - t331) * struc12PP + (-t24 * t12 * 
     &(-t577 + t206 + t426) + t3 * (t12 * (t18 * (t420 + t399 + t400 - t
     &199 + t39) - t209 - t211 + t231 - t424 - t504) + t41) - t331 - t12
     & * (t42 - t431 - t434 - t435 - t233 + t324 + t493) + t45) * struc1
     &3PP - t597 * struc14PP + t598 * struc15PP + t19 * struc17PP + t624
     & * struc21PP - t11 * struc22PP + t576 * struc24PP + t29 * struc25P
     &P + t32 * struc27PP + (t46 * t38 * t27 + t3 * (t12 * t40 + t41 * t
     &27) + t12 * t43 + t45 * t27) * struc16PP * t483)
      t3 = t7 * t197 * (-t114 * t143 * t456 * t256 + t24 * (t486 * t14 *
     & t606 * t12 + t8 * (t404 + t370)) + t3 * t10 + t46 * (t464 * (t26 
     &* (t456 * t457 * t191 - t456 * t37 * t382 - t639) - t328 * t456 * 
     &t606) + t549 * (-t9 * t22 + t24) + (t406 + t367 + t369) * t44 * s1
     &2 + t456 * ((t565 - t524) * s16 - t529 - t566) * t519 * propW16) +
     & t5) * struc35PP
      result = t174 * t2 + t6 * struc10PP + t1 * struc19PP + (-t126 * t1
     &74 + t166 * t207 * (t18 * (t54 + t399) - t209 - t234) + t169 * t56
     &8 * (t464 * (t12 * (s13 * t462 - t125 * t460) + t567 * t160) + (s1
     &6 * (t18 * (-t154 + t548) + t48 - t524) - t529 - t566) * t519 * pr
     &opW16) - t188 * t207 * t615 + t197 * t83 - t518) * struc20PP + t33
     & * struc30PP + t64 * struc31PP + t13 * struc32PP + (-t169 * t468 *
     & (t12 * (t460 * t60 + t462 * t53) + t26 * (t381 * t463 + t466)) - 
     &t166 * t471 * (-t234 + t469) - t174 * t74 + t197 * t57 + t188 * t4
     &70 * (t50 + t73) + t471 * (t163 * t233 - t231 * t89)) * struc33PP 
     &+ t4 * struc34PP + t105 * struc36PP + t3

           ampNonresonantLightFullPP = result
       end function ampNonresonantLightFullPP

       function ampNonresonantLightFullReC4MM()
           implicit none
           type(dd_complex) :: ampNonresonantLightFullReC4MM

           type(dd_complex) :: result

      result = ddreal(0)

           ampNonresonantLightFullReC4MM = result
       end function ampNonresonantLightFullReC4MM

       function ampNonresonantLightFullReC4MP()
           implicit none
           type(dd_complex) :: ampNonresonantLightFullReC4MP
           type(dd_complex) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           type(dd_complex) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           type(dd_complex) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           type(dd_complex) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           type(dd_complex) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           type(dd_complex) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           type(dd_complex) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           type(dd_complex) ::  t175,t176,t177,t178,t179,t18,t180,t181,t182,t183,t184,t185
           type(dd_complex) ::  t186,t187,t188,t189,t19,t190,t191,t192,t193,t194,t195,t196
           type(dd_complex) ::  t197,t198,t199,t2,t20,t200,t201,t202,t203,t204,t205,t206
           type(dd_complex) ::  t207,t208,t209,t21,t210,t211,t212,t213,t214,t215,t216,t217
           type(dd_complex) ::  t218,t219,t22,t220,t221,t222,t223,t224,t225,t226,t227,t228
           type(dd_complex) ::  t229,t23,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239
           type(dd_complex) ::  t24,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249,t25
           type(dd_complex) ::  t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t26,t260
           type(dd_complex) ::  t261,t262,t263,t264,t265,t266,t267,t268,t269,t27,t270,t271
           type(dd_complex) ::  t272,t273,t274,t275,t276,t277,t278,t279,t28,t280,t281,t282
           type(dd_complex) ::  t283,t284,t285,t286,t287,t288,t289,t29,t290,t291,t292,t293
           type(dd_complex) ::  t294,t295,t296,t297,t298,t299,t3,t30,t300,t301,t302,t303
           type(dd_complex) ::  t304,t305,t306,t307,t308,t309,t31,t310,t311,t312,t313,t314
           type(dd_complex) ::  t315,t316,t317,t318,t319,t32,t320,t321,t322,t323,t324,t325
           type(dd_complex) ::  t326,t327,t328,t329,t33,t330,t331,t332,t333,t334,t335,t336
           type(dd_complex) ::  t337,t338,t339,t34,t340,t341,t342,t343,t344,t345,t346,t347
           type(dd_complex) ::  t348,t349,t35,t350,t351,t352,t353,t354,t355,t356,t357,t358
           type(dd_complex) ::  t359,t36,t360,t361,t362,t363,t364,t365,t366,t367,t368,t369
           type(dd_complex) ::  t37,t370,t371,t372,t373,t374,t375,t376,t377,t378,t379,t38
           type(dd_complex) ::  t380,t381,t382,t383,t384,t385,t386,t387,t388,t389,t39,t390
           type(dd_complex) ::  t391,t392,t393,t394,t395,t396,t397,t398,t399,t4,t40,t400
           type(dd_complex) ::  t401,t402,t403,t404,t405,t406,t407,t408,t409,t41,t410,t411
           type(dd_complex) ::  t412,t413,t414,t415,t416,t417,t418,t419,t42,t420,t421,t422
           type(dd_complex) ::  t423,t424,t425,t426,t427,t428,t429,t43,t430,t431,t432,t433
           type(dd_complex) ::  t434,t435,t436,t437,t438,t439,t44,t440,t441,t442,t443,t444
           type(dd_complex) ::  t445,t446,t447,t448,t449,t45,t450,t451,t452,t453,t454,t455
           type(dd_complex) ::  t456,t457,t458,t459,t46,t460,t461,t462,t463,t464,t465,t466
           type(dd_complex) ::  t467,t468,t469,t47,t470,t471,t472,t473,t474,t475,t476,t477
           type(dd_complex) ::  t478,t479,t48,t480,t481,t482,t483,t484,t485,t486,t487,t488
           type(dd_complex) ::  t489,t49,t490,t491,t492,t493,t494,t495,t496,t497,t5,t50
           type(dd_complex) ::  t51,t52,t53,t54,t55,t56,t57,t58,t59,t6,t60,t61
           type(dd_complex) ::  t62,t63,t64,t65,t66,t67,t68,t69,t7,t70,t71,t72
           type(dd_complex) ::  t73,t74,t75,t76,t77,t78,t79,t8,t80,t81,t82,t83
           type(dd_complex) ::  t84,t85,t86,t87,t88,t89,t9,t90,t91,t92,t93,t94
           type(dd_complex) ::  t95,t96,t97,t98,t99

           type(dd_complex) :: result

      t1 = intHLs160000x0112D2eps0()
      t2 = intHLs16s25s26s34s56x1111D2eps0()
      t3 = intHs160000x0112D2eps0()
      t4 = intHs160s26s34s56x1021D2eps1()
      t5 = intHs160s26s34s56x1022D4eps0()
      t6 = intHs16s25s26s34s56x1112D4eps0()
      t7 = intHs16s25s26s34s56x1121D4eps0()
      t8 = intHs16s25s26s34s56x1211D4eps0()
      t9 = intHs16s25s26s34s56x1131D4eps0()
      t10 = ddreal(4)
      t11 = propZ25 * s25
      t12 = -t10 + t11
      t13 = intHs16s25s26s34s56x1221D4eps0()
      t14 = intHs16s25s26s34s56x1311D4eps0()
      t15 = intHs16s25s26s34s56x1311D4eps1()
      t16 = ddreal(2)
      t17 = -t16 + epinv
      t18 = ddreal(3)
      t19 = (gb)**2
      t20 = (gw)**2
      t21 = t18 * s25
      t22 = t21 * t20 * propZ25
      t23 = -t19 * t12 + t22
      t24 = ddreal(1) - epinv
      t25 = t19 * t12
      t26 = -t22 + t25
      t27 = intHLs16s25s26s34s56x1112D4eps0()
      t28 = t11 + t16
      t29 = t19 * t28
      t30 = t22 + t29
      t31 = intHs16s25s26s34s56x1310D4eps0()
      t32 = s16 + s26 - s34 + s56
      t33 = s12 + s25
      t34 = t10 * s26
      t35 = t16 * t33
      t36 = t35 + t34 + s16
      t37 = intHLs16s25s26s34s56x1211D2eps1()
      t38 = ddreal(6)
      t39 = ddreal(5)
      t40 = t18 * s26
      t41 = t38 * s12
      t42 = t39 * s25
      t43 = t40 + t42 + t41 + s34 - s56
      t44 = intHLs16s25s26s34s56x1121D2eps1()
      t45 = s25 + s26 + s56
      t46 = intHLs16s25s26s34s56x1112D2eps1()
      t47 = -s13 + s26 + s45 + s46 + s56
      t48 = t16 * t47 + s25
      t49 = intHs16s25s26s34s56x1222D6eps0()
      t50 = intHLs16s25s26s34s56x1121D2eps0()
      t51 = intHs16s25s26s34s56x1120D2eps0()
      t52 = intHs16s25s26s34s56x1210D2eps0()
      t53 = intHs16s25s26s34s56x1210D2eps1()
      t54 = intHLs160000x0111D0eps0()
      t55 = intHLs160000x0111D0eps1()
      t56 = intHLs160000x0112D2eps1()
      t57 = intHLs16s25s26s34s56x1112D4eps1()
      t58 = t16 * epinv
      t59 = -t18 + t58
      t60 = intHs16s25s26s34s56x1122D4eps0()
      t61 = t16 * s26
      t62 = s25 + t61
      t63 = intHs16s25s26s34s56x1212D4eps0()
      t64 = intHLs16s25s26s34s56x1222D6eps0()
      t65 = intHs16s25s26s34s56x1212D4eps1()
      t66 = t10 * s12
      t67 = intHLs16s25s26s34s56x1211D2eps0()
      t68 = intHLs16s25s26s34s56x1212D4eps1()
      t69 = s13 - s45 - s46
      t70 = t18 * s56
      t71 = t16 * t69
      t72 = -t70 + t66 + s25 + s26 - s34 + t71
      t73 = epinv * t72
      t74 = intHLs16s25s26s34s56x1212D4eps0()
      t75 = intHs16s25s26s34s56x1122D4eps1()
      t76 = t17 * s25
      t77 = intHLs16s25s26s34s56x1112D2eps0()
      t78 = intHLs16s25s26s34s56x1122D4eps1()
      t79 = t16 * s56
      t80 = -t69 + t79
      t81 = intHs16s25s26s34s56x1411D6eps0()
      t82 = intHs160s26s34s56x1031D4eps1()
      t83 = intHs16s25s26s34s56x1121D2eps0()
      t84 = intHs16s25s26s34s56x1131D4eps1()
      t85 = t76 - t66
      t86 = intHs16s25s26s34s56x1221D4eps1()
      t87 = t16 * s12
      t88 = t76 - t87
      t89 = intHs16s25s26s34s56x1141D6eps0()
      t90 = t87 + s25
      t91 = intHs16s25s26s34s56x1321D6eps0()
      t92 = t21 + t87
      t93 = intHs16s25s26s34s56x1231D6eps0()
      t94 = t21 + t66
      t95 = intHLs16s25s26s34s56x1132D6eps0()
      t96 = intHLs16s25s26s34s56x1312D6eps0()
      t97 = s25 + s16
      t98 = t16 * (s12 + s26)
      t99 = t98 + t97
      t100 = intHs16s25s26s34s56x1312D6eps0()
      t101 = t21 + t98
      t102 = intHs16s25s26s34s56x1113D4eps1()
      t103 = intHs16s25s26s34s56x1114D6eps0()
      t104 = intHs16s25s26s34s56x1213D6eps0()
      t105 = s12 + s25 + s26
      t106 = t105 * t16 + s34
      t107 = intHLs16s25s26s34s56x1113D4eps0()
      t108 = -t70 - s26 + t71
      t109 = intHLs16s25s26s34s56x1122D4eps0()
      t110 = intHLs16s25s26s34s56x1114D6eps0()
      t111 = s25 + s26 - s56
      t112 = intHLs16s25s26s34s56x1123D6eps0()
      t113 = t16 * s25
      t114 = t40 - s56 + t113
      t115 = intHLs16s25s26s34s56x1213D6eps0()
      t116 = t40 + t35 - s56
      t117 = intHs16s25s26s34s56x1132D6eps0()
      t118 = t66 + s26 + s34 - s56 + t113
      t119 = intHs16s25s26s34s56x1123D6eps0()
      t120 = t16 * (s12 + s34) + t111
      t121 = intHs16s25s26s34s56x1310D4eps1()
      t122 = intHs16s25s26s34s56x1110D2eps0()
      t123 = intHs160s26s34s56x1020D2eps0()
      t124 = intHLs16s25s26s34s56x1113D4eps1()
      t125 = intHs16s25s26s34s56x1220D4eps1()
      t126 = intHs16s25s26s34s56x1130D4eps1()
      t127 = intHs16s25s26s34s56x1321D6eps1()
      t128 = intHs16s25s26s34s56x1231D6eps1()
      t129 = intHLs16s25s26s34s56x1132D6eps1()
      t130 = intHLs16s25s26s34s56x1312D6eps1()
      t131 = intHs16s25s26s34s56x1312D6eps1()
      t132 = intHs16s25s26s34s56x1213D6eps1()
      t133 = intHLs16s25s26s34s56x1123D6eps1()
      t134 = intHLs16s25s26s34s56x1213D6eps1()
      t135 = intHs16s25s26s34s56x1132D6eps1()
      t136 = intHs16s25s26s34s56x1123D6eps1()
      t137 = intHs16s25s26s34s56x1110D2eps1()
      t138 = intHs16s25s26s34s56x1220D4eps0()
      t139 = intHs16s25s26s34s56x1130D4eps0()
      t140 = intHs16s25s26s34s56x1120D2eps1()
      t141 = t24 * s25
      t142 = intHs160s26s34s56x1020D2eps1()
      t143 = intHs160s26s34s56x1031D4eps0()
      t144 = ddreal(1) / s25
      t145 = t38 * t20 * propZ25
      t146 = t19 * (-propZ25 * t16 + ddreal(8) * t144) + t145
      t147 = intHs16s25s26s34s56x1211D2eps1()
      t148 = intHs160s26s34s56x1022D4eps1()
      t149 = intHs160000x0112D2eps1()
      t150 = intHs16s25s26s34s56x1121D2eps1()
      t151 = ddreal(1) + epinv
      t152 = t87 * t151
      t153 = intHs16s25s26s34s56x1121D4eps1()
      t154 = intHs16s25s26s34s56x1211D4eps1()
      t155 = intHs16s25s26s34s56x1112D4eps1()
      t156 = ddreal(1) / t32
      t157 = t138 * t33 + t139 * t90
      t158 = t125 * t33 + t126 * t90
      t159 = t62 * t129
      t160 = t99 * t130
      t161 = t101 * t131
      t162 = t106 * t132
      t163 = t118 * t135
      t164 = t120 * t136
      t165 = t6 * t144
      t166 = t2 * t28
      t167 = t166 * t144
      t168 = t15 * t17
      t169 = t27 * t30
      t170 = t24 * t128
      t171 = (t170 - t93) * t94
      t172 = t100 * t101
      t173 = t102 * s34
      t106 = t104 * t106
      t174 = t117 * t118
      t120 = t119 * t120
      t175 = t107 * t108
      t176 = -t109 * t80
      t177 = t62 * t95
      t178 = t96 * t99
      t179 = t84 * t85
      t180 = s12 * t83
      t181 = t24 * t121
      t182 = t142 * t24
      t183 = t17 * t137
      t184 = t153 + t154
      t185 = epinv * t149
      t186 = t24 * t148
      t187 = t150 * (t152 + s25)
      t188 = t7 * t12
      t189 = t48 * t77
      t190 = t59 * t155
      t191 = (t26 * (t144 * (-t182 + t123) + t51 + t52) + (t140 * (t141 
     &+ t87) + t122 + t183) * t144 * t23) * t156
      t192 = t144 * t30
      t193 = t24 * intHs16s25s26s34s56x1222D6eps1()
      t194 = t24 * intHLs16s25s26s34s56x1222D6eps1()
      t195 = t194 - t64
      t196 = t195 * t36
      t197 = (epinv * t37 + t67) * s34
      t198 = (epinv * t44 + t50) * t45
      t72 = t74 * t72
      t199 = t57 * t59
      t200 = (t60 + t63) * t62
      t201 = t75 * (epinv * t61 - s26 - s34 + s56 - t66 + t76)
      t202 = t53 * t24
      t203 = t202 * t156
      t204 = (t49 - t193) * t43
      t205 = t1 * t28
      t206 = t8 * t12
      t207 = t26 * (t144 * t204 - t203) + t143 * t146 + t23 * (t144 * (-
     &t186 + t187 + t185) + t147) + t191 + t144 * (-(t3 + t4 + t5) * t12
     & + t188) * t19 + t192 * (-t78 * (t16 * (epinv * t80 + s26) + s25) 
     &- t54 + t189) + t144 * (t17 * t184 + t190) * t26 + t192 * ((t46 * 
     &t48 - t55 + t56) * epinv + t68 * (t73 - s25 - t61) + t196 + t197 +
     & t198 + t72 - t199) + t144 * (t65 * ((-t18 + epinv) * s25 - t66 - 
     &t61 * t24) + t200 + t201) * t23 + t144 * (t206 + t205) * t19
      t208 = t24 * intHLs16s25s26s34s56x1114D6eps1()
      t209 = -t110 + t208
      t210 = t24 * intHs16s25s26s34s56x1141D6eps1()
      t211 = -t210 + t89
      t212 = t24 * intHs16s25s26s34s56x1114D6eps1()
      t213 = intHs16s25s26s34s56x1411D6eps1() * t24
      t214 = (-t212 + t103) * s34
      t215 = t211 * t90
      t216 = (t9 + t13 + t14) * s25
      t217 = t20 * propZ25
      t218 = intHL0s25s26s34s56x1110D2eps0()
      t219 = t45 - t71
      t220 = intHLs160000x0111D2eps0()
      t221 = t217 * s25
      t222 = t19 + t221
      t223 = s25 + s16 - s34
      t224 = t223 + t98
      t225 = t16 * (s12 - s45 - s46)
      t226 = t225 + s26 + s16 - s34 + s13 - s56
      t227 = intHLs16s25s26s34s56x1411D6eps0()
      t228 = intHs160s26s34s56x1011D2eps0()
      t229 = s26 + s16 - s34 - s56
      t230 = t16 * (s12 + s13 - s45 - s46)
      t231 = t230 + t229
      t232 = t16 * (s12 + s13 - s45) + s26 + s16 - s34 - s46 - s56
      t233 = t16 * s45
      t234 = s25 + s26 - s13 + s46 + s56 + t233
      t235 = intHL0s25s260s56x1012D2eps1()
      t236 = intHL0s25s260s56x1013D4eps1()
      t237 = intHL0s25s260s56x1010D0eps0()
      t238 = intHL0s25s260s56x1011D2eps0()
      t239 = intHL0s25s260s56x1020D2eps0()
      t240 = intHLs160000x0122D4eps0()
      t241 = intHs160s26s34s56x1011D2eps1()
      t242 = intHs160s26s34s56x1013D4eps0()
      t243 = intHL0s25s260s56x1011D2eps1()
      t244 = intHL0s25s26s34s56x1110D2eps1()
      t245 = intHL0s25s260s56x1010D0eps1()
      t246 = intHL0s25s260s56x1020D2eps1()
      t247 = intHL0s25s26s34s56x1120D2eps0()
      t248 = intHL0s25s26s34s56x1210D2eps0()
      t249 = intHL0s25s26s34s56x1310D4eps0()
      t250 = -s13 + s25 + s26 + s56
      t251 = s45 + s46
      t252 = t16 * t251
      t253 = t252 + t250
      t254 = intHs160s26s34s56x1012D2eps0()
      t255 = -t231
      t256 = intHs160s26s34s56x1021D2eps0()
      t257 = s12 + s13 + s16
      t258 = t16 * s13
      t259 = s12 - s26
      t260 = s13 - s45
      t261 = s12 + s13
      t262 = t261 * s13
      t263 = (-s34 + t257) * s25 + (s25 - s12 - s13) * s26 + (-s25 - s12
     & + s26 - t258 + s45) * s45 + (-t16 * t260 + s46 - t259) * s46 + (-
     &t259 - t71 + s56) * s56 + t262
      t264 = s25 - s12
      t265 = s25 - s13
      t266 = s26 - s13
      t267 = t16 * t266
      t268 = t16 * (s26 - s13 + s45)
      t269 = s25 - s12 + s16 - s34
      t270 = -s12 + s16 - s13
      t271 = s16 * s25
      t272 = t265 * s34
      t273 = (s56)**2
      t274 = t16 * t273 - t270 * s13 + (-t258 + t264) * s26 - (-t267 + s
     &12 - s16 + s34 - s45) * s45 - (-t223 - t268 - s46) * s46 + (-t18 *
     & t69 + t269 + t61) * s56 + t271 - t272
      t275 = t16 * s16
      t276 = s16 - s34 - s13
      t277 = t16 * t276 + t264 + t40
      t278 = s25 + s26
      t279 = (s34)**2
      t280 = (t277 + s45) * s45
      t281 = t261 * s25
      t260 = -t16 * (s16 * t261 - (s12 + s13 + s16 + s26) * s34 - t273) 
     &+ (-t18 * t261 - s26 - t275) * s26 + (t268 + t269 + s46) * s46 + (
     &-t10 * t260 - t18 * (s12 - s46) + t278) * s56 + t262 - t279 - (s16
     &)**2 + t280 - t281
      t269 = t18 * s13
      t282 = -t269 + t264
      t283 = t258 - s25
      t284 = -t10 * t69 + t16 * t223 - s12 + t40
      t285 = (t16 * (s25 + s16 - s34 - s13 + s45) + s12 + t40 + s46) * s
     &46
      t270 = -t270 * s25
      t286 = t283 * s34
      t287 = t18 * t273
      t288 = t258 * s16
      t280 = s26 * t282 + s56 * t284 + t262 - t270 + t280 + t285 + t286 
     &+ t287 - t288
      t289 = s16 + s26 + s56
      t290 = t18 * s34
      t291 = -t290 + t289
      t292 = t16 * s34
      t293 = -t292 + t289
      t294 = -t252 + t229
      t295 = s13 * t291 + t293 * t294 + t293 * t87
      t296 = intHs16s25s26s34s56x1112D2eps1()
      t297 = t16 * (s13 - s45 - s46 - s56)
      t298 = -t297 + s25
      t299 = intHL0s25s26s34s56x1220D4eps1()
      t300 = t45 + t292
      t301 = s25 + s26 + s34 + s56
      t252 = t252 + t45
      t302 = s13 * t300 - t252 * t301
      t303 = intHL0s25s260s56x1022D4eps1()
      t304 = s26 + s56
      t305 = t304 + t113
      t306 = s45 + s46 + s56
      t307 = t16 * t306
      t308 = s13 * t305
      t309 = t304 * t306
      t310 = (-t307 - s26 - s25) * s25 + t308 - t309
      t311 = intHLs160000x0121D2eps0()
      t312 = propW16 * s16
      t313 = t312 * t224
      t314 = t16 * (-t313 + s25)
      t315 = t10 * s45
      t316 = s46 + s56
      t317 = t313 * t38
      t318 = t19 * (t28 * (-t18 * t316 - s26 - t113 + t269 - t315) + t31
     &7) + t22 * (t18 * (s13 - s46 - s56) - s26 - t314 - t315)
      t319 = s13 - s45 - s56
      t320 = t18 * s46
      t314 = t10 * t319 - t314 - t320
      t321 = t10 * s13
      t322 = t28 * (t321 - t10 * (s45 + s56) - t113 - t320) + t317
      t323 = t18 * epinv
      t324 = s25 + s56
      t325 = t16 * t304
      t307 = t307 + s25
      t326 = (t21 + t325) * t307 - t258 * (t16 * t324 + s26)
      t327 = t16 * t32
      t328 = (t10 * t32 + t87) * s13 - t307 * (-s25 + t327)
      t329 = intHLs160000x0121D2eps1()
      t330 = intHLs16s25s26s34s56x1211D4eps1()
      t331 = s25 - s26
      t332 = s16 - s34
      t333 = s13 - s56
      t334 = t38 * s45
      t335 = t39 * s46
      t336 = t16 * t332
      t337 = t18 * t333
      t338 = -t336 + t335 + t334 - t66 - t337 + t331
      t339 = ddreal(7)
      t340 = t18 * t332
      t341 = t10 * t333
      t342 = ddreal(8) * s45
      t343 = t339 * s46
      t344 = t38 * t251
      t264 = t16 * t264
      t345 = t264 + t344 - t337 + s26 - s16 + s34
      t259 = t16 * t259
      t346 = ddreal(8) * t251
      t347 = intHLs16s25s26s34s56x1121D4eps1()
      t334 = t316 * t39 + s26 + t21 - t321 + t334
      t348 = t10 * s25
      t349 = t39 * s13
      t350 = s16 + s26 - s34
      t351 = s56 * t39 + t264 - t321 + t344 - t350
      t352 = t339 * s56
      t335 = -t319 * t38 + t21 + t335
      t353 = s26 + s45 + s56
      t354 = s26 * s46
      t355 = (s25)**2
      t356 = s46 * s56
      t357 = t349 * s25
      t344 = t16 * (-s13 * t304 + s26 * s45 + s56 * t353) + t18 * (t356 
     &+ t355) + s25 * (t344 + t352 + s26) + t354 - t357
      t358 = t18 * s16
      t359 = t16 * t261
      t360 = -s16 + s25 + s34
      t361 = s12 + s16
      t362 = t16 * t361
      t363 = -t269 + s26 + t113
      t364 = (s26)**2
      t365 = t251 * t360
      t366 = s12 * s13
      t367 = t363 * s34
      t368 = s13 * s16
      t369 = (-s12 + s13 - s16 - s25 - s26) * s26
      t370 = (-t358 - t359 + s25) * s25
      t371 = (t21 - t321 + t61) * s34
      t372 = t16 * (-(t336 - s25 + s12 + s26) * s56 + t369 + t366)
      t373 = t10 * (t368 + t365)
      t374 = s26 - s34
      t375 = (s46)**2
      t376 = (-s12 - s16 - s25 - s26) * s26
      t377 = (-s13 + s25 + s26) * s34
      t378 = t257 * s25
      t257 = t16 * (s13 * t257 + (s25 - s16 + s34 - s13) * s45 - t375 + 
     &t376 + t377 - t378) - (t16 * (s16 + s45) + t374) * s46 - (t16 * (s
     &12 + s26 + s16 - s34 + s13) + t320) * s56
      t379 = -t21 + t321
      t380 = s25 - s26 - s16 + s34
      t381 = s26 + s16 - s34 - s13 + s45 + s46
      t382 = t251 * t380
      t383 = (t18 * t265 + s16 + s26 + t87) * s26
      t276 = (t16 * (s12 + s45 + s46) + t18 * t276 + s25 + t34) * s56
      t384 = t16 * (t366 + t382 - t378)
      t385 = t18 * (t368 - t273)
      t386 = t87 * s13
      t387 = intHLs16s25s26s34s56x1221D4eps1()
      t388 = s13 + s16
      t389 = t87 + t388
      t390 = -s13 + s16
      t391 = s25 + s26 + s34
      t392 = s25 + s45 + s46
      t393 = t10 * t392
      t394 = t39 * s34
      t395 = t40 - s13 + t113
      t396 = s12 - s45
      t397 = t16 * t396
      t398 = s25 - s16 + s34 - s13 + s46 - t397
      t399 = (t87 + s25 + s16 + s13 + s26) * s26
      t400 = t391 * s46
      t401 = (s16 + t359) * s25
      t402 = t16 * (-s45 * t391 + t366)
      t403 = (-t87 - t390 - s26) * s26
      t404 = (-t389 + s25) * s25
      t266 = -t113 - t266
      t405 = (t349 + t66 + t275) * s25
      t233 = (s25 - t233) * s46
      t406 = (s13)**2
      t407 = (-s25 + s12 - s45 + s46 + t258 - s56) * s56
      t408 = (s12 + s13 + s16 + s25 + s26) * s26
      t409 = t278 * s34
      t410 = t16 * t381 + s25
      t411 = s26 + s16
      t267 = (-s25 - t267) * s34 - t16 * ((-t16 * t411 - s12 + s13 - s25
     & + s34 - s45 - s46 - s56) * s56 - t251 * (t374 + t275) - t271 + t3
     &66 + t369) - t321 * s16
      t271 = t251 * t350
      t369 = s13 * t361
      t412 = t369 - t271 - t273
      t413 = t390 * s25
      t414 = -t283 * s26
      t415 = -t349 + s16 + t113
      t416 = s13 * t339 - s26 - t348
      t417 = s26 + s45 + s46
      t418 = t339 * s34
      t419 = t10 * t417
      t420 = t18 * t390
      t261 = t10 * t261
      t421 = t16 * (s26 + s16 - s13 + s45 + s46)
      t422 = t251 * (t411 - t292)
      t279 = -t16 * (t368 - t422 + t279 - t273) + (t388 + s26) * s25 + (
     &t40 + t261 + t275) * s34 + (-t290 + s25 + t421) * s56
      t390 = t113 * t390
      t423 = t16 * t283
      t424 = s25 - s34
      t271 = -t16 * (t368 - t271 - t273) + s25 * t388 + (t411 + t359 - s
     &34) * s34 + (t424 + t421) * s56 + t414
      t359 = t396 * s25
      t396 = (-t353 + s13) * s13
      t229 = t229 * s25
      t354 = t16 * (-s45 * s46 - t359 - t375 + t396) - t18 * (s13 * s25 
     &+ t356) - t354 - t229
      t356 = s25 * t251
      t421 = s26 - s56
      t425 = s25 * s26
      t426 = t10 * t356
      t229 = t16 * (-s46 * t353 - t359 - t375 + t396) - t229 - t21 * s13
      t353 = s26 + s34
      t359 = t16 * (s34 * t306 + t355) + t426 + s13 * (-t21 + t336) - s1
     &6 * t307 + t353 * s25 + t21 * s56
      t396 = intHLs16s25s26s34s56x1141D6eps0()
      t427 = s13 * (s25 + t79)
      t428 = t45 * t307
      t429 = t427 - t428
      t430 = t392 * s16
      t431 = s12 * t45
      t432 = intHLs16s25s26s34s56x1131D4eps1()
      t433 = s26 * s34
      t434 = s25 + t327
      t435 = t32 - t113
      t436 = s13 * (-t21 + t327) - t307 * t435
      t331 = t16 * (t366 + t355) - t287 - (t358 + s13) * s25 - (-t349 + 
     &t97 + s26) * s26 - s34 * t416 + (-s16 * t39 + s25 + t269 + t418 - 
     &t419) * s56 + t251 * (-s16 * t38 + t10 * t331 + ddreal(8) * s34) +
     & t349 * s16
      t437 = t62 * s34
      t438 = intHLs16s25s26s34s56x1321D6eps0()
      t439 = s25 + s16 - s13
      t440 = t10 * s34
      t264 = -t16 * (t251 * (t292 + t278) + t366) + (t87 + s16) * s25 - 
     &(-t439 - t87 - s26) * s26 - (-t269 + t264 - s16 + s34) * s34 - (-t
     &225 + t440 + s25 - s16 - s13 + s56) * s56 - t368
      t417 = t16 * (-t251 * (-t290 + t411) + t369) + (t269 - t97 - s26) 
     &* s26 + (-t16 * t417 - s56 - t439 + t440) * s56 + t371 - t413
      t439 = intHLs16s25s26s34s56x1131D4eps0()
      t441 = t16 * t278
      t442 = intHLs16s25s26s34s56x1221D4eps0()
      t443 = intHLs16s25s26s34s56x1231D6eps0()
      t444 = t16 * (t251 * (t441 + s34) + t366) + t287 + (t16 * t265 + s
     &26) * s34 + (t18 * (s34 - s13) - s16 - t259 + t393) * s56 + t368 +
     & t403 + t404
      t445 = t339 * s25
      t446 = t10 * t251
      t265 = t10 * (t353 + t113) * t251 - t16 * (-s16 * t69 - t355 - t36
     &6) + t287 - (t358 + t261) * s25 + (t10 * t265 + s26) * s34 + (-t18
     & * t388 - t259 + t394 + t445 + t446) * s56 - t399
      t353 = intHs160000x0211D2eps0()
      t388 = intHL0s25s260s56x1031D4eps0()
      t447 = -s16 + s25
      t448 = -t447 - t290
      t449 = t10 * s34 * t251 + t16 * (s34 + s45 + s46) * s25 + s13 * t4
     &48 - s16 * t252 + s25 * t45 + t433 + t290 * s56
      t450 = intHL0s25s260s56x1022D4eps0()
      t451 = -t16 * (s12 + s26 + s13 - s45 - s46 - s56) - t332
      t452 = intHL0s25s260s56x1021D2eps1()
      t453 = t24 * s26
      t454 = epinv * s56
      t455 = intHLs16s25s26s34s56x1231D6eps1()
      t456 = t24 * s34
      t457 = -t456 + s26 + s16 + s56
      t458 = t58 * s34
      t459 = t458 + t32
      t460 = intHs160000x0211D2eps1()
      t461 = -t10 + t323
      t462 = intHL0s25s260s56x1031D4eps1()
      t463 = intHLs16s25s26s34s56x1311D4eps1()
      t464 = intHL0s25s26s34s56x1220D4eps0()
      t465 = intHL0s25s26s34s56x1120D2eps1()
      t466 = epinv * s34
      t467 = -t466 + s25 + s26 + s56
      t468 = -t16 * t319 + s25 + s46
      t411 = t10 * t369 - t287 - t390 + (-t415 - s26) * s26 + (t21 - t34
     &9 + s26) * s34 + (t394 - t419 - t420 - t113) * s56 + t251 * (s34 *
     & t38 - t10 * t411)
      t469 = intHLs16s25s26s34s56x1321D6eps1()
      t470 = intHL0s25s26s34s56x1130D4eps1()
      t471 = intHs16s25s26s34s56x1113D4eps0()
      t472 = intHLs160000x0211D2eps1()
      t473 = t16 * (s12 + s13 - s56)
      t315 = t320 - t473 + t315 + t360
      t320 = -t473 + t446 + t360
      t473 = intHs160000x0121D2eps1()
      t380 = t446 - t87 - t337 + t380
      t292 = -t447 - t292
      t389 = t16 * t292 * t251 - (t87 + s13 + s16) * s16 + s25 * t389 + 
     &s26 * t447 - (t16 * (s25 - s12 - s16) - s26 - t269 + s34) * s34 + 
     &s56 * t448
      t272 = t16 * (t368 + t272) - (s26 + t269 + t362 - s25) * s25 + (-t
     &336 + t348) * s45 - (t336 - t348 + s26) * s46 + (t21 - t336 + s46)
     & * s56
      t448 = intHL0s25s260s56x1012D2eps0()
      t474 = intHL0s25s260s56x1013D4eps0()
      t475 = intHL0s25s260s56x1021D2eps0()
      t476 = intHs160000x0111D0eps0()
      t477 = intHL0s25s26s34s56x1130D4eps0()
      t478 = intHLs160000x0211D2eps0()
      t479 = intHs160000x0121D2eps0()
      t480 = intHLs16s25s26s34s56x1211D4eps0()
      t481 = intHLs16s25s26s34s56x1121D4eps0()
      t482 = intHs16s25s26s34s56x1111D2eps0()
      t483 = t223 + t61
      t484 = propW16 * t224
      t485 = t484 * t45 * t32
      t486 = s12 + s16 + s26 - s34
      t487 = intHs16s25s26s34s56x1211D2eps0()
      t488 = intHL0s25s26s34s56x1210D2eps1()
      t489 = intHL0s25s26s34s56x1310D4eps1()
      t490 = intHs160s26s34s56x1012D2eps1()
      t491 = intHs160s26s34s56x1013D4eps1()
      t492 = ddreal(1) / t45
      t493 = t51 + t52
      t494 = ddreal(1) / (ecossin)**2
      t495 = s34 * t488
      t496 = epinv * t55
      t497 = t24 * t492
      t41 = -t330 * (epinv * t338 - s25 + t340 + t341 - t342 - t343 + t4
     &1 + t61) + t46 * (epinv * t229 + s25 * t253) - t57 * (epinv * t335
     & + ddreal(8) * t319 - t343 - t348) - t68 * (epinv * t257 + t370 + 
     &t371 + t372 + t373) - t78 * (t10 * (t251 * (s26 + t113) + t273 + t
     &355) - t357 + epinv * (-t16 * t266 * s45 - t16 * (-t409 + t408 + t
     &407 + t375 - t406) + t233 - t405) + (t21 - t258) * s26 + (ddreal(9
     &) * s25 + t10 * (s26 - s13 + s45 + s46)) * s56) - t347 * (epinv * 
     &t334 - t316 * t339 - s26 - t342 - t348 + t349) - t387 * (t10 * t39
     &1 * t251 + t287 + t386 + epinv * (s34 * t395 - (-t398 - s56) * s56
     & - t399 + t400 - t401 - t402) + (t21 - t321 + s26) * s34 + (-t269 
     &- t259 + t394 + t393 - s16) * s56 + t368 + t403 + t404) - t496 * t
     &468 + t497 * (t219 * t495 + t299 * t302)
      t316 = t123 * t156
      t319 = intHs160000x0111D0eps1() * epinv
      t286 = t150 * (epinv * t267 + s56 * t410 - t16 * t412 + t286 + t41
     &3 + t414) - t75 * (-t10 * (t369 - t422) + t287 + epinv * t279 + (t
     &415 + s26) * s26 + s34 * t416 + (t419 - t418 + t420 + t113) * s56 
     &+ t390) + t456 * t255 * t490 * t156
      t321 = epinv * t329
      t342 = t494 * t144 * propW34
      t343 = epinv * t472
      t357 = s34 * intHs16s25s26s34s56x1112D2eps0()
      t369 = t154 * t461
      t391 = t144 * t26
      t149 = t30 * (t144 * (t194 * t265 - t257 * t74 + t315 * t343) + t4
     &92 * (t144 * (t219 * (t218 - t237 + t238 + t239) - t302 * t464 - t
     &452 * (-t16 * (s13 * (-t141 - s26 - s56) + t309) - s25 * (t141 - t
     &16 * (-s45 * t24 - s46 * t24) + t70 - t454 + t453)) - t465 * (s13 
     &* (-t458 + t45) - t252 * t467)) + t24 * t219 * t235) + t144 * (t37
     & * (epinv * t468 + s25 - t297) + t468 * t67) * s34) + t144 * ((-t1
     &49 * t320 - t380 * t473) * epinv + t156 * (-t122 * t232 - t228 * t
     &231) - t482 * (s13 + s26 - s56) - t357 * t298) * t23 + t391 * (-t1
     &56 * (t140 * (epinv * t263 - s26 * t282 - t283 * s34 - (t277 + s45
     &) * s45 - s56 * t284 - t262 + t270 - t285 - t287 + t288) + t295 * 
     &t5 + t4 * (s13 * t459 + t294 * t457 + t457 * t87)) - t369 * t298)
      t252 = t494 * propW34
      t250 = -t13 * (-t363 * s34 + t276 + t383 - t384 - t385) - t279 * t
     &60 - t3 * t320 - t63 * (t16 * (s12 * t250 + s26 * t304 - s34 * t39
     &2 - t356 + t430) + t18 * (s56 * t332 + t425) + s13 * (-t340 + s25 
     &- t61) - s25 * s56 + s26 * t332) - t380 * t479
      t257 = t30 * (t448 + t475) * t492
      t229 = t252 * (t144 * (t23 * t250 + t26 * (-t253 * t476 - t267 * t
     &83 - t331 * t49 - t345 * t6 - t351 * t7) + t30 * (t2 * (s25 + s26 
     &+ s16 - s34 - s13 - s46 + s56 + t87) - t27 * t335 - t334 * t481 - 
     &t338 * t480 - t44 * (-epinv * (-t16 * (s13 * (s25 - s16 + s56) - s
     &25 * s45 + (s25 + s16 - s34 - s45) * s26 + (-s25 + s16 - s45) * s5
     &6 - t273 + t364 + t430 + t431) + s25 * s34 + s46 * t45) + t427 - t
     &428) - t50 * (-t16 * ((s25 + s26 - s16) * s45 + t273 + t368 + t376
     & - t378) + (t275 - t278) * s46 + (-t16 * (s25 - s12 - s16 - s13 + 
     &s45) - s46) * s56 - t437) - t54 * t468 - t64 * t265 + t77 * t229 +
     & t442 * (-s34 * t395 + (-t398 - s56) * s56 + t399 - t400 + t401 + 
     &t402) + t315 * t478)) - t257 * t219)
      t250 = (t156)**2
      t262 = (t243 + t244) * t17
      t265 = t492 * t219
      t267 = t151 * s34
      t270 = t487 * t23
      t69 = t144 * (propW34 * (t23 * (t17 * (-t137 * t232 - t231 * t241)
     & - s34 * t255 * (t254 + t256)) + (t148 * t295 + t263 * t53) * t26 
     &* t24) - epinv2 * (t19 * (-t16 * propW34 * (t16 * (-t251 * (t40 + 
     &t332 + t113) + t368) - t287 + (t261 + s16) * s25 + (s13 * t38 + s2
     &6 + t66 + t97) * s26 - (t258 + t278) * s34 + (t38 * t69 - t21 + t3
     &32 - t61 + t66) * s56 + t11 * (-s56 * t486 + t365 + t368 + t376 + 
     &t377 - t378)) + t485 * t38 * propW34) - t145 * propW34 * s25 * (s1
     &3 * t332 + (t258 + s12) * s26 + (-s25 + s12 - s26 + t71 - s56) * s
     &56 - t251 * t483 + t281 - t485)) * t492) * t156
      t69 = propW34 * (-t270 * t253 + t147 * t26 * (epinv * t253 - s13 -
     & s16 - s26 + s34 + s56 - t225)) + t144 * (-t311 * propW34 * t318 -
     & t1 * propW34 * (t19 * t322 + t22 * t314) - t56 * (t19 * propW34 *
     & (epinv * t322 - t317) + t221 * propW34 * (t314 * t323 - t317)) - 
     &t267 * propW34 * t23 * t296 * t298) - t182 * t391 * propW34 * t260
     & * t250 + t69 + t192 * propW34 * (t109 * (t16 * (s45 * t266 + t375
     & - t406 + t407 + t408 - t409) + t405 - t233) + t265 * ((-t247 - t2
     &48) * s34 + t24 * (t245 - t246) + t262))
      t97 = t463 * s34
      t145 = epinv * t460
      t225 = s34 * t471
      t232 = s34 * intHLs16s25s26s34s56x1311D4eps0()
      t233 = t24 * t236
      t251 = t252 * (t144 * (-t23 * (t271 * t9 + t298 * t353) + t26 * (-
     &t100 * t359 - t104 * t449 - t117 * t411 - t143 * t226 + t24 * (t12
     &7 * t436 + t128 * t328 + t226 * t82) - t173 * (epinv * t451 + s26 
     &+ t113 - t337 + t446)) + t30 * (-t107 * t354 - t115 * t272 + t24 *
     & (t129 * t326 + t133 * t344) + t253 * t477 + t298 * t388 - t432 * 
     &(-epinv * (-t16 * (s45 * t45 - t431 + t433) - (-s16 + s34 + s46) *
     & s25 - (-s25 - s16 + s46) * s26 - (s25 - s16 + s46) * s56 - t273 +
     & t308 + t364) - t16 * t429) + t96 * t389 + t439 * (-s46 * t278 - (
     &s25 - s16 - s13 + s46 - t397 + s56) * s56 + t399 + t401 - t437 - t
     &441 * s45) - t443 * t444 - t232 * t234 - t450 * t310 * t492)) + t2
     &65 * t30 * (-t233 + t474))
      t255 = t252 * (t26 * (-t14 * t298 - t144 * (t119 * t417 + t436 * t
     &91) - t15 * (epinv * t298 - s16 + s25 - s26 + s34 - t337 + t446 - 
     &t87)) + t192 * (s34 * t249 * t265 + t264 * t438) + t144 * (t24 * (
     &t121 * t263 + t125 * t274 + t126 * t280) - s34 * t231 * t242) * t1
     &56 * t23)
      t219 = t255 + t251 + t342 * (t23 * (t156 * (t231 * t456 * t491 - t
     &138 * t274 - t139 * t280 - t263 * t31) - t84 * (s56 * t16 * t410 +
     & epinv * t271 - s26 * t423 + s34 * t423 - t10 * t412 + t390) - t14
     &5 * t298) + t26 * (t24 * (t131 * t359 + t132 * t449 + t135 * t411 
     &+ t136 * t417) - t328 * t93 - t225 * t451) + t30 * (-t112 * t344 -
     & t124 * (-s25 * t18 * t333 + epinv * t354 - s46 * t421 + t16 * t35
     &5 + t425 + t426) + t24 * (-t130 * t389 + t134 * t272 - t253 * t470
     & - t264 * t469 - t298 * t462 + t444 * t455 + t492 * (-s34 * t219 *
     & t489 + t303 * t310)) - t326 * t95 - t97 * (epinv * t234 - s16 + s
     &25 - s26 + s34 - t337 + t446 - t87)))
      t231 = t24 * intHLs160000x0113D4eps1() - intHLs160000x0113D4eps0()
      t234 = t24 * intHLs160000x0122D4eps1()
      t251 = intHLs160000x0111D2eps1() * t17
      t255 = (t234 - t240) * s16 - t220 - t251
      t110 = -t208 + t110
      t89 = -t210 + t89
      t208 = -t212 + t103
      t210 = t24 * intHLs16s25s26s34s56x1411D6eps1()
      t258 = t24 * intHLs16s25s26s34s56x1141D6eps1()
      t226 = t252 * (t144 * (t26 * ((s13 * t434 - t307 * t32 + t386) * t
     &89 + t208 * t253 * s34) + t30 * (t429 * (t258 - t396) + t110 * (-t
     &16 * s25 * (-t306 + s13) - s46 * t421 + t355) + (t210 - t227) * t2
     &26 * s34)) + t484 * t144 * t255 * t222 + t26 * (-t213 + t81) * t29
     &8)
      t261 = -ddreal(1) / ddreal(3)
      t41 = t10 * t342 * t313 * t222 * t231 + t16 * t226 + t261 * (t494 
     &* t69 + t229 + t252 * t149 + t342 * (t23 * t286 + t26 * (-t153 * (
     &epinv * t351 + s16 + s26 - s34 - t21 - t346 + t349 - t352 + t87) -
     & t155 * (epinv * t345 + s16 - s34 - t21 + t259 + t341 - t346) + t1
     &56 * (t260 * t316 - t263 * t493) - t65 * (t372 + t373 + epinv * (t
     &16 * (t366 + t365 - t364) - (t362 + s13) * s25 - (t16 * (s12 - s13
     &) + s16 + t21) * s26 + (-t340 - t98 + s25) * s56 + t367 + t269 * s
     &16) + t370 + t371) - t86 * (t10 * (t368 + t382 - t273) + epinv * (
     &t385 + t384 - t383 + t367 - t276) - (t358 + t87 + s13 - s25) * s25
     & + s26 * t379 - t379 * s34 + (-t10 * t381 + s25) * s56 + t386) - t
     &319 * t253 + t193 * t331) + t30 * t41 - t321 * t318)) - ddreal(2) 
     &/ ddreal(3) * t219 + t342 * t8 * t26 * t298
      t69 = -s16 + s25 - s26 + s34 - s56
      t149 = s25 + s34
      t219 = epinv * t32
      t226 = epinv * s25
      t229 = t57 * t461
      t253 = t12 * t479
      t259 = epinv * t473
      t260 = t125 * t69
      t261 = t260 * t156
      t263 = t303 * t305
      t264 = epinv * t4
      t265 = t24 * t490
      t148 = t156 * (t23 * (t148 * t456 - t228) + s34 * (t265 - t264) * 
     &t26) + t23 * (t13 * t435 + t86 * (epinv * t435 + s16 + s26 - s34 +
     & s56)) + t30 * (t301 * t64 - t32 * t442 - t387 * (t219 - s34) + t4
     &92 * (s34 * (epinv * t465 + t247) + t17 * (-t243 - t244) + t24 * (
     &(t299 - t488) * s34 - t245 + t246 + t263) - t218 - t238 - t239) + 
     &t78 * (epinv * t149 + s25 + s26 + s56)) - t49 * t26 * t293
      t243 = t140 * (t226 + s26 + s16 - s34 + s56)
      t244 = t450 * t305
      t266 = t452 * (t141 + s26 + s56)
      t137 = (t241 + t137) * t17
      t241 = t138 * t69
      t269 = t144 * (-(t254 + t256) * s34 + t182 - t243)
      t271 = epinv * t150
      t272 = t65 * t151 * t23
      t273 = t192 * ((-s34 * t46 - t44 * t447) * epinv + t492 * ((-t464 
     &+ t248) * s34 - t244 - t266) + t54 - t194 * t301 - t267 * t37)
      t274 = -t122 * t23 - t123 * t26
      t276 = s34 * t67
      t277 = t166 * t19
      t279 = t30 * (t448 + t475) * t492
      t148 = t144 * (t156 * t274 + t30 * (t481 - t276 + t496) - t277) + 
     &t144 * (t26 * (t193 * t293 - t267 * t296 - t271 * t447) - t272 * t
     &360) + t156 * (t26 * (t269 - t51) - epinv2 * (t19 * ((propZ25 * t3
     &60 - t10) * s25 - t304 * t38 - t336) - t22 * (t223 + t325)) * t144
     & * t492 + t144 * (-s34 * t5 - t137 - t241) * t23) + t144 * t148 + 
     &t144 * (t19 * ((t482 + t7) * t12 - t28 * t311 - t253) + t26 * (-t4
     &47 * t83 - t259 - t357)) + t30 * (t144 * (t109 * t149 + t237 * t49
     &2 - t447 * t50 + t68 * (s25 - s16 + t267) + t229 - t321) + t24 * (
     &-t144 * t347 - t235 * t492) + t144 * (t74 - t77) * s34) + t144 * (
     &t24 * (t261 + t153) + t60 * t292 - t360 * t63 - t75 * (-epinv * t2
     &92 + s34) - t369) * t23 + t273 + t26 * t156 * (-t52 + t202) + t279
      t149 = -t213 + t81
      t223 = t149 * t23
      t267 = t209 * t30
      t273 = t24 * t491
      t279 = t24 * t126
      t280 = t121 * t156
      t281 = t31 * t156
      t69 = t19 * (t12 * t14 + t144 * (-t12 * (t139 + t353) - t205)) + t
     &23 * (t144 * (-epinv * t69 * t84 - s34 * t156 * t242) - t15 * t151
     & - t281) + t24 * (t23 * (t144 * (s34 * t132 + t127 * t69 + t131 * 
     &t360) + t280) + t391 * t128 * t32 - t192 * (s34 * t130 + t129 * t4
     &5 + t133 * t305 + t134 * t360)) + t30 * ((-t144 * t56 + t432) * ep
     &inv + t492 * (-s34 * t144 * t249 - t474)) + t144 * (t26 * (-t32 * 
     &t93 - t145 + t279) - t29 * t388) + t30 * (t144 * ((t489 * t497 + t
     &107 + t463 + t96) * s34 + t112 * t305 + t115 * t360 + t124 * (t466
     & + s25) + t24 * t462 + t45 * t95) + t439 + t497 * t236) + t144 * (
     &(-t102 * t151 + t156 * t273 - t104 - t471) * s34 - t100 * t360 + t
     &69 * (-t9 - t91)) * t23
      t145 = ddreal(9) * t217
      t236 = t38 * propW16
      t282 = t312 * t16
      t283 = ddreal(1) + t282
      t284 = s25 + s26 - s34 + s56
      t285 = t312 * t38
      t286 = t285 + t11 + t16
      t287 = t460 + t473
      t288 = t7 + t8
      t292 = t1 * t286
      t293 = s34 * t296
      t294 = epinv2 * (t19 * (t236 * t32 - t10 + t11) + t22 * (propW16 *
     & t327 - ddreal(1))) * t156
      t286 = t56 * ((t19 * t286 + t22 * t283) * epinv - t285 * t222)
      t295 = t46 * epinv
      t298 = -t267 + ((t240 - t234 - t311 - t321) * s16 + t220 + t251) *
     & t144 * t222 * propW16 + t217 * (-t216 + t139)
      t242 = -t273 + t242
      t273 = t24 * t133
      t302 = t273 - t112
      t303 = t24 * t134
      t304 = t303 - t115
      t306 = t129 * t24
      t307 = t306 - t95
      t308 = t130 * t24 - t96
      t309 = t156 * t23
      t305 = t144 * (t30 * (t308 * s34 + t107 * t424 - t124 * (t141 + t4
     &66) + t305 * t302 + t360 * t304 + t45 * t307 - t27) - t279 * t26 +
     & s34 * (t156 * t242 + t102) * t23) + t25 * (t139 * t144 - t13 - t1
     &4 - t9) + (t84 + t86 + t15) * t23 * epinv + t309 * (-t181 + t31)
      t283 = t1 * t283
      t310 = t312 * t222
      t313 = t310 * t144 * t231
      t195 = t16 * t305 + t18 * t217 * (t283 + t2 + t7 + t8 - t353 - t47
     &9) - t38 * t298 + ddreal(12) * t313 - t156 * (t144 * t274 + t202 *
     & t26) + t144 * (t26 * t271 * t284 + t25 * t353 + t272 * t45 + t286
     &) - t156 * (t26 * (t144 * ((-t264 + t265 - t254 - t256) * s34 + t1
     &82 - t243) - t51 - t52) + t144 * ((-t5 + t186) * s34 - t137 - t241
     &) * t23) - t144 * (t19 * (t12 * t288 - t253 - t292) + t23 * (-t156
     & * t228 + t24 * (t261 + t153 + t154) + t45 * (-t60 - t63) - t75 * 
     &(epinv * t45 + s34)) + t26 * (-epinv * t287 - t284 * t83 - t293) +
     & t30 * (-s34 * t77 - t109 * t424 + t289 * t74 + t68 * (epinv * t28
     &9 + s25 + s26 + s56) + t78 * (t141 + t466 + s26 + s56)) - t294) - 
     &t192 * (-s34 * t295 - t195 * t301 + t54) + t167 * t19 - t192 * (-t
     &197 - t198 + t496 + t199)
      t243 = t16 + epinv
      t253 = -s13 + s25 + s26 + s45 + s56
      t254 = s12 + s26 + s16 - s34 + s56
      t256 = t16 * t254
      t261 = t10 * t361
      t254 = t10 * t254
      t264 = t254 + s25
      t271 = t312 * t10
      t272 = ddreal(1) + t271
      t274 = t16 * (s12 + s26 + s16) - t290
      t284 = -t16 * t486 + s25
      t289 = -t256 + s25
      t290 = t40 + t362 - t440 + t324
      t298 = -t16 * (s16 - s34 + s13 - s45 - s56) + s25 - t66
      t305 = t339 * t374 + t361 * t38 + s25 + t70
      t314 = t16 * (s26 - s13 + s45 + s56) + t21
      t315 = t10 * s56
      t268 = t21 + t315 + t268
      t317 = t39 * s26
      t312 = ddreal(12) * t312
      t318 = t312 + t11 + t16
      t135 = t24 * t135
      t320 = t24 * t82
      t322 = t45 * t439
      t323 = t432 * t45 + t97
      t324 = t144 * t24
      t325 = t15 * t243
      t326 = t323 * epinv
      t327 = t309 * t122 * t144
      t197 = t192 * (-t197 - t198 + t496 + t177 + t178)
      t117 = t26 * (t144 * (t136 * t24 * t290 - t316) + t203) + t156 * (
     &t26 * (t269 - t51 - t52) + t144 * (-t241 - t137) * t23) + t23 * (t
     &324 * (t127 * t289 + t131 * t284) - t325) + t144 * (t19 * (-t12 * 
     &(t143 + t479) - t292) + t23 * (-t100 * t284 + t156 * (t24 * t260 -
     & t228) - t361 * t83 - t84 * ((-t16 + epinv) * s25 - t254) - t86 * 
     &(t226 - t256) - t91 * t289) + t26 * ((-t460 - t473) * epinv + (t15
     &6 * t265 - t296) * s34 - t104 * t274 - t119 * t290 - t93 * t264 + 
     &(t135 - t117) * (t374 * t39 + t113 + t261 + t70) + t320) + t30 * (
     &t124 * (t253 * t58 + s25) + t232 + t322) - t294) - t286 * t144 + t
     &192 * (t24 * (-t159 - t160) + t54 + t326) + t25 * (-t144 * t353 + 
     &t13 + t14 + t9) + t391 * t24 * (t128 * t264 + t132 * t274) - t327 
     &+ t197
      t119 = (t330 + t347) * t24
      t128 = t301 * t442
      t136 = t301 * t387
      t36 = t144 * (t26 * (t156 * (t291 * t5 + t4 * t459) - t193 * t305)
     & + t30 * ((t314 * t46 - t136) * epinv + t36 * (t194 - t64) - t480 
     &- t481) + t277 + (t60 + t63) * t62 * t23) + t146 * t6 + t26 * (t14
     &4 * (-t156 * t186 * t291 + t305 * t49 + t185) + t147) + t321 * (t1
     &9 * t318 + t22 * t272) * t144 + t192 * (-t109 * t268 - t298 * t74 
     &+ t314 * t77 - t68 * (epinv * t298 + s25 + t61) - t78 * (epinv * t
     &268 + s25 + t61) - t229 + t119 + t343 - t128) + t144 * (t150 * (t1
     &6 * (t151 * t361 + s26 - s34 + s56) + s25) + t184 * t461 + t65 * (
     &s25 * t151 - t10 * (s12 + s16 - s34) - t16 * t453) + t75 * (t16 * 
     &(epinv * s26 - s25) + t226 - t261 + t418 - t70 - t317) + t190) * t
     &23 + (t12 * t3 + t28 * t478 + t311 * t318) * t144 * t19
      t137 = t7 + t8
      t190 = t25 * t137
      t36 = -t10 * (t144 * (t26 * (t279 - t139) + t30 * (-t105 * t304 + 
     &t107 * t253 - t278 * t302)) + t309 * (-s34 * t144 * t242 + t181 - 
     &t31)) - t117 * t16 + t18 * (t144 * (-t169 - t190) + t217 * (t272 *
     & t311 + t2 - t3 + t478)) + ddreal(24) * t313 + t38 * ((t102 + t103
     & - t212) * t23 * t144 * s34 + t217 * (t283 - t479 + t216 - t143 - 
     &t353) + t391 * t211 * (t256 + s25) + t267 + t223) + t36 - ddreal(1
     &2) * t144 * t222 * propW16 * ((-t234 + t240) * s16 + t220 + t251) 
     &+ ddreal(9) * t217 * t137
      t102 = s25 - s26 + s56
      t103 = s16 + s56
      t117 = t21 + t317 + t103 + t87
      t197 = t40 + t35 + s16 + s34 + s56
      t79 = s25 - s13 + s45 + s46 + t79
      t198 = t70 - t71 + t278
      t211 = epinv * t48
      t212 = t151 * s12
      t220 = t16 * (t212 + s25 + s26)
      t228 = t21 + s26 - s16 + s34 - s56 + t87
      t229 = t40 + t348 + t87 + s56
      t55 = t19 * (-(t4 + t5) * t12 + t205) + t23 * (t60 * t62 + t201) +
     & t30 * ((s34 * t37 - t55 + t56) * epinv + t276)
      t28 = t144 * t55 + t23 * (t144 * (t259 - t186 + t187) + t147 * t15
     &1) + t191 + t144 * (t26 * (t13 * (t32 + t87) + t153 * t59 + t155 *
     & t461 + t224 * t63 + t49 * t228 + t65 * (epinv * t483 + t220) + t8
     &6 * (t219 + t21 + t152) + t319) + t30 * (t102 * t50 + t17 * t347 +
     & t24 * (t57 + t299) + t330 * t59 - t387 * (-t24 * t87 - s16 - s56 
     &- t21 + t226 - t317) - t44 * (s26 * t243 + t141 - t454) - t442 * t
     &90 + t46 * (t211 - s25) + t64 * t229 + t68 * (epinv * (-t70 + t374
     & + t71) + t220) + t74 * (-t70 + t230 + t374) - t78 * (t16 * (epinv
     & * t79 - s26) - t21))) + t146 * (t3 + t143) + t270 - t391 * t193 *
     & t228 + t192 * (-t194 * t229 + t189 - t27 - t343 - t464 - t465 - t
     &478 - t54) + t144 * (t12 * t476 - (t479 + t482) * t12 + t28 * t481
     &) * t19 - t203 * t26
      t55 = t144 * (t26 * (t215 + t214) + t30 * ((-t210 + t227) * s34 + 
     &t209 * t421 + t62 * (-t258 + t396))) + t217 * (s25 * t9 + t2 + t48
     &0 - t7)
      t71 = t165 * t25 + t217 * (-t476 + t479 + t481 + t482 + t1 + t4 + 
     &t5)
      t15 = t10 * t144 * (t30 * t278 * (-t306 + t95) + (-t170 + t93) * t
     &33 * t26) - t16 * (t19 * (t12 * t9 + t144 * (-t188 - t166)) + t23 
     &* (t144 * (-t156 * t157 - t179 - t180 - t185) - t281) + t24 * (t23
     & * (t144 * t82 + t280) + t192 * (t102 * (t133 + t134) + t117 * t45
     &5 + t130 * t360 + t197 * t469) + t391 * (t132 * t360 + t163 + t164
     &)) + t26 * (t24 * (t127 + t131) - t100 - t15 - t91 - t104 * t360 *
     & t144) + t30 * (t144 * (t107 * t198 + t109 * t79 + t124 * (epinv *
     & t198 - s25 + s26 - s56) - t24 * t470 - t360 * t96 + t432 * (t76 -
     & t34) - t438 * t197 - t443 * t117 - (t112 + t115) * t102 + t477 - 
     &t97) + t439) + t324 * t309 * t158 - t29 * t480 * t144 - t391 * (t1
     &73 + t174 + t120)) + t18 * t71 + t38 * t55 + t28 - t145 * t6
      t28 = -s13 + s25 + s45 + s46 + s56
      t29 = -s12 - s13 + s25 + s45 + s46 + s56
      t33 = s12 - s13 + s25 + s45 + s46 + s56
      t34 = s12 + s13 + s16 + s26 - s34 - s45 - s46
      t55 = t10 * t374 + t103 * t16 + s25
      t71 = t16 * t45
      t70 = t70 + t348 - s26
      t42 = t42 + t315 + t61
      t76 = s25 * t38 + t10 * (s34 + s56) - t358 - t87
      t79 = t45 + t275
      t87 = -t16 * t350 + s25
      t97 = -t358 + t440 - t98 + s25
      t98 = t19 * (-t285 + t11 + t16) + t22 * (ddreal(1) - t282)
      t103 = ddreal(1) + t58
      t104 = t19 * (t312 - t11 - t16) + t22 * (-ddreal(1) + t271)
      t117 = t91 * t92
      t11 = t342 * (t23 * (t13 * t32 + t434 * t9 + t84 * (epinv * t434 -
     & t113 - t66)) + t26 * (t86 * (-t219 + t35) + t93 * t94 + t106 + t3
     &19 + t117 - t118 * t135) + t30 * (t112 * t70 + t124 * (epinv * t10
     &8 - s26 + s56 + t21) + t492 * (-t301 * t464 - t218 + t237 - t238 -
     & t239) + t95 * t42 + t96 * t97 + (-t303 + t115) * (t16 * (s25 - s1
     &2 - s16 + s34) + s56 - t40) + t175) + epinv2 * (t19 * (t236 * t45 
     &+ t11 + t16) + t22 * (propW16 * t71 + ddreal(1))) * t492)
      t22 = t92 * t127
      t35 = t252 * (t144 * (t23 * (t353 + t479) + t26 * (-t170 * t94 + t
     &476) + t30 * (t2 + t27 - t478) - t311 * t98) + t257 + t325 * t26)
      t11 = t35 + t342 * (t23 * (t156 * (t138 * t28 + t139 * t33 + t29 *
     & t31) + t143 + t293 * epinv) + t26 * (-t24 * (t22 + t161 + t162 + 
     &t164) + t173 * t103) + t30 * (t24 * (-t129 * t42 - t130 * t97 - t1
     &33 * t70 + t492 * (-t245 + t246)) + t492 * ((t247 + t248) * s34 - 
     &t465 * t467 - t244 - t266) - t262 * t492)) + t11 + t252 * (t23 * (
     &t144 * (t83 * (s12 - s16 - s26 + s34 - s56) + t357) + t487) + t24 
     &* (t144 * t23 * (-t156 * (t121 * t29 + t125 * t28 + t126 * t33) - 
     &t82) + (t144 * (t299 * t301 + t263 - t495) - t235) * t492 * t30) +
     & t26 * (t144 * (t120 + t172 + t174) + t14) + t144 * (t23 * t287 - 
     &t30 * t472 - t329 * t98) * epinv + t192 * (t59 * (t330 + t347) + t
     &176))
      t28 = t252 * (t144 * (-t1 * t104 - t23 * (t156 * t183 + t4 + t5) -
     & t56 * (epinv * t104 - ddreal(12) * t310)) - t147 * t103 * t23 + t
     &391 * (t156 * (t140 * (t24 * t28 + t212) - t29 * t493 - t316 * t34
     &) - t63 * t87 - t65 * (epinv * t87 + t21 + t61 + t66) + t193 * t43
     &) + t192 * (-(t17 * t37 + t67) * s34 - t46 * (t211 - t113) - t68 *
     & (-t10 * t332 + t21 - t61 + t73) + t54 + t194 * t76))
      t28 = t28 + t342 * (t23 * (-t122 * t156 - t150 * (-t16 * (t219 - t
     &212) + s25) - t55 * t60 - t75 * (epinv * t55 - s26 - s34 + s56 - t
     &113 - t66) + t186) + t26 * (-t155 * (epinv * t10 - t39) - t43 * t4
     &9 + t461 * (-t153 - t154) + t24 * t156 * (t142 * t156 * t34 + t29 
     &* t53)) + t30 * (-t44 * (epinv * t79 - t71) - t50 * t79 - t64 * t7
     &6 - t78 * (t16 * (-epinv * t80 + s26) + t315 + t445) - t189 - t199
     & + t496 - t72))
      t29 = t297 + s25
      t33 = epinv * t62
      t34 = t297 + t66 + s25
      t35 = t143 - t6
      t22 = t23 * (t180 - t320 - t106 - t174 - t120) + t30 * (t232 + t32
     &2) + t23 * ((t168 + t9 + t13 + t14) * s25 + t101 * (t131 * t24 - t
     &100) + t24 * (t162 + t163 + t164 + t22) + t84 * t85 + t86 * t88 - 
     &t117 + t171) + t30 * (-t124 * (t47 * t58 + s25) + t307 * t62 + t30
     &8 * t99 + t326) + t173 * t26 - t25 * t35
      t25 = t252 * ((-t16 * t69 - t18 * (t144 * (t19 * t206 + t169) + t2
     &17 * (t479 - t482 - t2 - t7 - t311)) + t38 * (t267 + t223 + t217 *
     & (s25 * t14 + t1 - t139 - t353 + t388)) - t148 + t145 * t8) * stru
     &c2MP + (struc3MP + struc8MP) * t195 + t192 * (-t16 * (epinv * t124
     & + t107) + (t46 - t68 - t78) * epinv - t109 - t74 + t77) * struc61
     &MP)
      t17 = struc9MP * (t10 * t30 * (t105 * (t303 - t115) - t107 * t47 +
     & t278 * (t273 - t112)) + t16 * t22 - t18 * t221 * t137 + t38 * (t2
     &3 * (-s25 * t149 - t214 - t215) + t267 * s25 + t221 * t35) - t23 *
     & (-s25 * t147 + t59 * t155 + t17 * t184 + t75 * (t66 + t113 - t33 
     &+ s26 + s34 - s56) - t185 + t186 - t187 - t200) - t30 * ((t329 + t
     &472) * epinv - t109 * t29 + t57 * t17 + t119 - t128) + t190 - t23 
     &* (t65 * (-t33 + t21 + t66 + t61) + t204 - t3 - t4 - t5) - t30 * (
     &t48 * (-t295 - t77) + t68 * (-epinv * t34 + s25 + t61) - t74 * t34
     & + t78 * (-epinv * t29 + s25 + t61) - t196 - t2 + t27 + t311 + t47
     &8 - t480 - t481 - t136 * epinv)) - struc60MP * t26 * (t32 * (-t140
     & * t24 - t16 * (-t125 * t24 + t138 + t139 - t181 - t279 + t31) - t
     &202 + t51 + t52) - t123 + t182) * t250
      t1 = t252 * (t144 * t17 + (-t16 * (t23 * (-t144 * (t86 * t88 + t17
     &9 + t180) + t181 * t156) + t19 * (t12 * (-t165 + t9 + t13 + t14) -
     & t167) + t23 * (t156 * (-t144 * t157 - t31) - t168) + t169 * t144 
     &+ t144 * (t23 * (t156 * t158 + t82) + t26 * (t162 + t163 + t164 + 
     &t161) + t30 * (-t114 * t133 - t116 * t134 - t159 - t160)) * t24 + 
     &t144 * (t92 * (t127 * t24 - t91) + t171 - t172 - t173 - t106 - t17
     &4 - t120) * t26 + (t112 * t114 + t115 * t116 + t124 * (-epinv * t1
     &08 + s25 + s26 - s56) - t175 - t176 + t177 + t178) * t144 * t30) +
     & t18 * t217 * (t1 + t3 + t4 + t5 - t7 - t8) + t38 * (t26 * (t144 *
     & t214 + t144 * t215 - t213 + t81) + t192 * t209 * t111 + (t216 + t
     &2 - t6) * propZ25 * t20) + t207) * struc10MP + t36 * struc5MP + t1
     &5 * struc6MP)
      result = ddreal(4) / ddreal(3) * t25 + t41 * struc1MP + (t10 * t25
     &2 * (t26 * (s34 * t144 * t208 + t144 * t89 * t90 - t213 + t81) + t
     &192 * t110 * t102) + t16 * t342 * (t23 * (t185 + t3) + t26 * t288)
     & + ddreal(4) / ddreal(3) * t11 + ddreal(16) * t310 * t342 * t231 +
     & ddreal(8) / ddreal(3) * t252 * (t30 * (t144 * (t24 * (-t300 * t46
     &9 + t462 + t470) + t300 * t438 + t136 - t388 - t477 + t480 + t481)
     & + t492 * ((t24 * t489 - t249) * t144 * s34 - t474 + t233) + t144 
     &* (-t24 * t455 + t443) * (s34 + t71)) + t391 * (t225 + t6)) + ddre
     &al(16) / ddreal(3) * t252 * t192 * t323 + ddreal(8) * t342 * (t30 
     &* ((-t210 + t227) * s34 + t45 * (-t258 + t396)) + t255 * t222 * pr
     &opW16) - ddreal(2) / ddreal(3) * t28) * struc7MP + ddreal(2) / ddr
     &eal(3) * t1

           ampNonresonantLightFullReC4MP = result
       end function ampNonresonantLightFullReC4MP

       function ampNonresonantLightFullReC4PM()
           implicit none
           type(dd_complex) :: ampNonresonantLightFullReC4PM
           type(dd_complex) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           type(dd_complex) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           type(dd_complex) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           type(dd_complex) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           type(dd_complex) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           type(dd_complex) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           type(dd_complex) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           type(dd_complex) ::  t175,t176,t177,t178,t179,t18,t180,t181,t182,t183,t184,t185
           type(dd_complex) ::  t186,t187,t188,t189,t19,t190,t191,t192,t193,t194,t195,t196
           type(dd_complex) ::  t197,t198,t199,t2,t20,t200,t201,t202,t203,t204,t205,t206
           type(dd_complex) ::  t207,t208,t209,t21,t210,t211,t212,t213,t214,t215,t216,t217
           type(dd_complex) ::  t218,t219,t22,t220,t221,t222,t223,t224,t225,t226,t227,t228
           type(dd_complex) ::  t229,t23,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239
           type(dd_complex) ::  t24,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249,t25
           type(dd_complex) ::  t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t26,t260
           type(dd_complex) ::  t261,t262,t263,t264,t265,t266,t267,t268,t269,t27,t270,t271
           type(dd_complex) ::  t272,t273,t274,t275,t276,t277,t278,t279,t28,t280,t281,t282
           type(dd_complex) ::  t283,t284,t285,t286,t287,t288,t289,t29,t290,t291,t292,t293
           type(dd_complex) ::  t294,t295,t296,t297,t298,t299,t3,t30,t300,t301,t302,t303
           type(dd_complex) ::  t304,t305,t306,t307,t308,t309,t31,t310,t311,t312,t313,t314
           type(dd_complex) ::  t315,t316,t317,t318,t319,t32,t320,t321,t322,t323,t324,t325
           type(dd_complex) ::  t326,t327,t328,t329,t33,t330,t331,t332,t333,t334,t335,t336
           type(dd_complex) ::  t337,t338,t339,t34,t340,t341,t342,t343,t344,t345,t346,t347
           type(dd_complex) ::  t348,t349,t35,t350,t351,t352,t353,t354,t355,t356,t357,t358
           type(dd_complex) ::  t359,t36,t360,t361,t362,t363,t364,t365,t366,t367,t368,t369
           type(dd_complex) ::  t37,t370,t371,t372,t373,t374,t375,t376,t377,t378,t379,t38
           type(dd_complex) ::  t380,t381,t382,t383,t384,t385,t386,t387,t388,t389,t39,t390
           type(dd_complex) ::  t391,t392,t393,t394,t395,t396,t397,t398,t399,t4,t40,t400
           type(dd_complex) ::  t401,t402,t403,t404,t405,t406,t407,t408,t409,t41,t410,t411
           type(dd_complex) ::  t412,t413,t414,t415,t416,t417,t418,t419,t42,t420,t421,t422
           type(dd_complex) ::  t423,t424,t425,t426,t427,t428,t429,t43,t430,t431,t432,t433
           type(dd_complex) ::  t434,t435,t436,t437,t438,t439,t44,t440,t441,t442,t443,t444
           type(dd_complex) ::  t445,t446,t447,t448,t449,t45,t450,t451,t452,t453,t454,t455
           type(dd_complex) ::  t456,t457,t458,t459,t46,t460,t461,t462,t463,t464,t465,t466
           type(dd_complex) ::  t467,t468,t469,t47,t470,t471,t472,t473,t474,t475,t476,t477
           type(dd_complex) ::  t478,t479,t48,t480,t481,t482,t483,t484,t485,t486,t487,t488
           type(dd_complex) ::  t489,t49,t490,t491,t492,t493,t494,t495,t496,t497,t498,t499
           type(dd_complex) ::  t5,t50,t500,t501,t502,t503,t51,t52,t53,t54,t55,t56
           type(dd_complex) ::  t57,t58,t59,t6,t60,t61,t62,t63,t64,t65,t66,t67
           type(dd_complex) ::  t68,t69,t7,t70,t71,t72,t73,t74,t75,t76,t77,t78
           type(dd_complex) ::  t79,t8,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89
           type(dd_complex) ::  t9,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99

           type(dd_complex) :: result

      t1 = intHLs160000x0112D2eps0()
      t2 = intHLs16s25s26s34s56x1111D2eps0()
      t3 = intHs160000x0112D2eps0()
      t4 = intHs160s26s34s56x1021D2eps1()
      t5 = intHs160s26s34s56x1022D4eps0()
      t6 = intHs16s25s26s34s56x1112D4eps0()
      t7 = intHs16s25s26s34s56x1121D4eps0()
      t8 = intHs16s25s26s34s56x1211D4eps0()
      t9 = intHs16s25s26s34s56x1131D4eps0()
      t10 = ddreal(4)
      t11 = propZ25 * s25
      t12 = t11 - t10
      t13 = intHs16s25s26s34s56x1221D4eps0()
      t14 = intHs16s25s26s34s56x1311D4eps0()
      t15 = intHs16s25s26s34s56x1311D4eps1()
      t16 = ddreal(2)
      t17 = -t16 + epinv
      t18 = ddreal(3)
      t19 = (gb)**2
      t20 = (gw)**2
      t21 = t18 * s25
      t22 = t21 * t20 * propZ25
      t23 = -t19 * t12
      t24 = t22 + t23
      t25 = ddreal(1) - epinv
      t26 = t19 * t12
      t27 = -t22 + t26
      t28 = intHLs16s25s26s34s56x1112D4eps0()
      t29 = t11 + t16
      t30 = t19 * t29
      t31 = t22 + t30
      t32 = intHs16s25s26s34s56x1310D4eps0()
      t33 = s16 + s26 - s34 + s56
      t34 = intHLs16s25s26s34s56x1222D6eps1()
      t35 = s12 + s25
      t36 = t16 * t35
      t37 = t10 * s26
      t38 = t37 + s16 + t36
      t39 = intHLs16s25s26s34s56x1211D2eps1()
      t40 = intHs16s25s26s34s56x1222D6eps1()
      t41 = ddreal(6)
      t42 = ddreal(5)
      t43 = t18 * s26
      t44 = t42 * s25
      t45 = t41 * s12
      t46 = t43 + s34 - s56 + t45 + t44
      t47 = intHLs16s25s26s34s56x1121D2eps1()
      t48 = s25 + s26 + s56
      t49 = intHLs16s25s26s34s56x1112D2eps1()
      t50 = -s13 + s26 + s45 + s46 + s56
      t51 = t16 * t50 + s25
      t52 = intHs16s25s26s34s56x1220D4eps1()
      t53 = intHs16s25s26s34s56x1130D4eps1()
      t54 = t16 * s12
      t55 = t54 + s25
      t56 = intHs16s25s26s34s56x1411D6eps0()
      t57 = intHLs16s25s26s34s56x1213D6eps0()
      t58 = t43 - s56 + t36
      t59 = intHs16s25s26s34s56x1132D6eps0()
      t60 = t16 * s25
      t61 = t10 * s12
      t62 = t60 + s26 + s34 - s56 + t61
      t63 = intHs16s25s26s34s56x1123D6eps0()
      t64 = s25 + s26 - s56
      t65 = t16 * (s12 + s34) + t64
      t66 = intHs16s25s26s34s56x1310D4eps1()
      t67 = intHs16s25s26s34s56x1110D2eps0()
      t68 = intHs160s26s34s56x1020D2eps0()
      t69 = intHLs16s25s26s34s56x1113D4eps1()
      t70 = s13 - s45 - s46
      t71 = t18 * s56
      t72 = t16 * t70
      t73 = t71 - t72 + s26
      t74 = intHs16s25s26s34s56x1321D6eps1()
      t75 = t54 + t21
      t76 = intHs16s25s26s34s56x1231D6eps1()
      t77 = t61 + t21
      t78 = intHLs16s25s26s34s56x1132D6eps1()
      t79 = t16 * s26
      t80 = t79 + s25
      t81 = intHLs16s25s26s34s56x1312D6eps1()
      t82 = s25 + s16
      t83 = t16 * (s12 + s26)
      t84 = t82 + t83
      t85 = intHs16s25s26s34s56x1312D6eps1()
      t86 = t83 + t21
      t87 = intHs16s25s26s34s56x1213D6eps1()
      t88 = s12 + s25 + s26
      t89 = t16 * t88 + s34
      t90 = intHLs16s25s26s34s56x1123D6eps1()
      t91 = t43 + t60 - s56
      t92 = intHLs16s25s26s34s56x1213D6eps1()
      t93 = intHs16s25s26s34s56x1132D6eps1()
      t94 = intHs16s25s26s34s56x1123D6eps1()
      t95 = intHs16s25s26s34s56x1110D2eps1()
      t96 = intHs16s25s26s34s56x1220D4eps0()
      t97 = intHs16s25s26s34s56x1130D4eps0()
      t98 = intHs16s25s26s34s56x1120D2eps1()
      t99 = t25 * s25
      t100 = intHs160s26s34s56x1020D2eps1()
      t101 = intHs16s25s26s34s56x1122D4eps1()
      t102 = t17 * s25
      t103 = t16 * epinv
      t104 = intHLs16s25s26s34s56x1112D2eps0()
      t105 = intHLs16s25s26s34s56x1122D4eps1()
      t106 = t16 * s56
      t107 = t106 - t70
      t108 = intHs16s25s26s34s56x1222D6eps0()
      t109 = intHLs16s25s26s34s56x1121D2eps0()
      t110 = intHs16s25s26s34s56x1120D2eps0()
      t111 = intHs16s25s26s34s56x1210D2eps0()
      t112 = intHs16s25s26s34s56x1210D2eps1()
      t113 = intHLs160000x0111D0eps0()
      t114 = intHLs160000x0111D0eps1()
      t115 = intHLs160000x0112D2eps1()
      t116 = intHLs16s25s26s34s56x1112D4eps1()
      t117 = t103 - t18
      t118 = intHs16s25s26s34s56x1122D4eps0()
      t119 = intHs16s25s26s34s56x1212D4eps0()
      t120 = intHLs16s25s26s34s56x1222D6eps0()
      t121 = intHs16s25s26s34s56x1212D4eps1()
      t122 = intHLs16s25s26s34s56x1211D2eps0()
      t123 = intHLs16s25s26s34s56x1212D4eps1()
      t124 = s25 + s26 - s34
      t125 = -t71 + t124 + t61 + t72
      t126 = intHLs16s25s26s34s56x1212D4eps0()
      t127 = intHs160000x0112D2eps1()
      t128 = intHs16s25s26s34s56x1121D2eps1()
      t129 = ddreal(1) + epinv
      t130 = t54 * t129
      t131 = intHs16s25s26s34s56x1121D4eps1()
      t132 = intHs16s25s26s34s56x1211D4eps1()
      t133 = intHs16s25s26s34s56x1112D4eps1()
      t134 = intHs160s26s34s56x1031D4eps0()
      t135 = ddreal(1) / s25
      t136 = t41 * t20 * propZ25
      t137 = t19 * (-propZ25 * t16 + ddreal(8) * t135) + t136
      t138 = intHs16s25s26s34s56x1211D2eps1()
      t139 = intHs160s26s34s56x1031D4eps1()
      t140 = intHs16s25s26s34s56x1121D2eps0()
      t141 = intHs16s25s26s34s56x1131D4eps1()
      t142 = t102 - t61
      t143 = intHs16s25s26s34s56x1221D4eps1()
      t144 = -t54 + t102
      t145 = intHs16s25s26s34s56x1141D6eps0()
      t146 = intHs16s25s26s34s56x1321D6eps0()
      t147 = intHs16s25s26s34s56x1231D6eps0()
      t148 = intHLs16s25s26s34s56x1132D6eps0()
      t149 = intHLs16s25s26s34s56x1312D6eps0()
      t150 = intHs16s25s26s34s56x1312D6eps0()
      t151 = intHs16s25s26s34s56x1113D4eps1()
      t152 = intHs16s25s26s34s56x1114D6eps0()
      t153 = intHs16s25s26s34s56x1213D6eps0()
      t154 = intHLs16s25s26s34s56x1113D4eps0()
      t155 = intHLs16s25s26s34s56x1122D4eps0()
      t156 = intHLs16s25s26s34s56x1114D6eps0()
      t157 = intHLs16s25s26s34s56x1123D6eps0()
      t158 = ddreal(1) / t33
      t159 = t17 * t95
      t160 = t131 + t132
      t161 = t25 * t100
      t162 = t117 * t133
      t163 = t25 * t40
      t164 = t49 * t51
      t165 = t39 * s34
      t166 = t47 * t48
      t167 = t25 * t34
      t168 = t1 * t29
      t169 = t128 * (t130 + s25)
      t170 = epinv * t127
      t171 = t25 * intHs160s26s34s56x1022D4eps1()
      t172 = (t98 * (t54 + t99) + t67 + t159) * t158
      t173 = t7 + t8
      t174 = t122 * s34
      t175 = t104 * t51
      t176 = t48 * t109
      t177 = t116 * t117
      t178 = t38 * t120
      t179 = (-t114 + t115) * epinv
      t180 = (t118 + t119) * t80
      t181 = t101 * (t103 * s26 - s26 - s34 + s56 + t102 - t61)
      t182 = t112 * t25
      t183 = t25 * intHLs16s25s26s34s56x1114D6eps1()
      t184 = t183 - t156
      t185 = t25 * intHs16s25s26s34s56x1141D6eps1()
      t186 = -t185 + t145
      t187 = t25 * intHs16s25s26s34s56x1114D6eps1()
      t188 = intHs16s25s26s34s56x1411D6eps1() * t25
      t189 = (-t187 + t152) * s34
      t190 = t186 * t55
      t191 = (t9 + t13 + t14) * s25
      t192 = t184 * t31
      t193 = t25 * t53
      t194 = -t193 + t97
      t195 = t25 * t52
      t196 = t65 * t94
      t197 = t74 * t75
      t198 = t76 * t77
      t199 = t85 * t86
      t200 = t87 * t89
      t201 = t62 * t93
      t202 = (t96 - t195) * t35
      t203 = t194 * t55
      t204 = t15 * t17
      t205 = t78 * t80
      t206 = t81 * t84
      t207 = t6 * t135
      t208 = t135 * t31
      t209 = t135 * t25
      t210 = t80 * t148
      t84 = t84 * t149
      t86 = t86 * t150
      t89 = t89 * t153
      t211 = s34 * t151
      t212 = t59 * t62
      t213 = t63 * t65
      t214 = t75 * t146
      t215 = t77 * t147
      t216 = s12 * t140
      t217 = t25 * t139
      t218 = t141 * t142
      t219 = t2 * t29
      t220 = t219 * t19
      t221 = t25 * t66
      t58 = t135 * (t24 * (t143 * t144 + t216 - t217 + t218) + t220) + t
     &24 * (t158 * (t135 * (t202 + t203) + t32) + t204) + t208 * (t25 * 
     &(t58 * t92 + t90 * t91 + t205 + t206) - t28) + t26 * (-t9 - t13 - 
     &t14 + t207) + t209 * (-t198 - t199 - t200 - t196 - t197 - t201) * 
     &t27 + (t215 + t211 + t212 + t213 + t214 + t86 + t89) * t135 * t27 
     &+ t208 * (-t57 * t58 - t69 * (epinv * t73 + s25 + s26 - s56) - t91
     & * t157 - t154 * t73 - t155 * t107 - t210 - t84) - t221 * t24 * t1
     &58
      t73 = t20 * propZ25
      t20 = t16 * t58 - t18 * t73 * (-t1 - t3 - t4 - t5 + t7 + t8) + t41
     & * (t27 * (t189 * t135 + t190 * t135 - t188 + t56) + t192 * t135 *
     & t64 + (t191 + t2 - t6) * propZ25 * t20) + t27 * (t46 * t108 * t13
     &5 - t182 * t158) + t134 * t137 + t135 * (t31 * ((t165 + t166 + t16
     &4) * epinv - t113 + t167 * t38) + t168 * t19) + t24 * (t135 * (t17
     &2 - t171 + t170 + t169) + t138) + t27 * (t135 * (-t163 * t46 + t16
     &2) + t158 * (t135 * (-t161 + t68) + t110 + t111) + t135 * t160 * t
     &17) + t135 * (-t105 * (t16 * (epinv * t107 + s26) + s25) + t123 * 
     &(epinv * t125 - s25 - t79) + t126 * t125 + t174 + t175 + t176 - t1
     &77 - t178 + t179) * t31 + t135 * (t121 * ((-t18 + epinv) * s25 - t
     &61 - t79 * t25) + t180 + t181) * t24 + t135 * (t12 * t173 - (t3 + 
     &t4 + t5) * t12) * t19
      t58 = intHL0s25s26s34s56x1110D2eps0()
      t64 = t48 - t72
      t91 = intHLs160000x0111D2eps0()
      t107 = t73 * s25
      t125 = t107 + t19
      t222 = s25 + s16 - s34
      t223 = t222 + t83
      t224 = s26 + s16 - s34 + s13 - s56
      t225 = t16 * (s12 - s45 - s46)
      t226 = t224 + t225
      t227 = intHLs16s25s26s34s56x1411D6eps0()
      t228 = intHs160s26s34s56x1011D2eps0()
      t229 = s26 + s16 - s34 - s56
      t230 = t16 * (s12 + s13 - s45 - s46)
      t231 = t229 + t230
      t232 = t16 * (s12 + s13 - s45) + s26 + s16 - s34 - s46 - s56
      t233 = s45 * t16 - s13 + s25 + s26 + s46 + s56
      t234 = intHL0s25s260s56x1010D0eps0()
      t235 = intHL0s25s260s56x1011D2eps0()
      t236 = intHL0s25s260s56x1020D2eps0()
      t237 = intHLs160000x0122D4eps0()
      t238 = intHs160s26s34s56x1011D2eps1()
      t239 = intHs160s26s34s56x1013D4eps0()
      t240 = intHL0s25s260s56x1011D2eps1()
      t241 = intHL0s25s26s34s56x1110D2eps1()
      t242 = intHL0s25s260s56x1010D0eps1()
      t243 = intHL0s25s260s56x1020D2eps1()
      t244 = intHL0s25s26s34s56x1120D2eps0()
      t245 = intHL0s25s26s34s56x1210D2eps0()
      t246 = intHL0s25s26s34s56x1310D4eps0()
      t247 = -s13 + s25 + s26 + s56
      t248 = s45 + s46
      t249 = t16 * t248
      t250 = t247 + t249
      t251 = intHs160s26s34s56x1012D2eps0()
      t252 = -t231
      t253 = intHs160s26s34s56x1021D2eps0()
      t254 = s12 + s13 + s16
      t255 = t16 * s13
      t256 = s12 - s26
      t257 = s13 - s45
      t258 = t16 * t257
      t259 = s12 + s13
      t260 = t259 * s13
      t261 = (-s34 + t254) * s25 + (s25 - s12 - s13) * s26 + (-t255 - s2
     &5 - s12 + s26 + s45) * s45 + (-t256 - t258 + s46) * s46 + (-t256 -
     & t72 + s56) * s56 + t260
      t262 = s25 - s12
      t263 = s25 - s13
      t264 = t16 * (s26 - s13)
      t265 = -s13 + s26 + s45
      t266 = t16 * t265
      t267 = s25 - s12 + s16 - s34
      t268 = s12 + s13 - s16
      t269 = t263 * s34
      t270 = s16 * s25
      t271 = (s56)**2
      t272 = t16 * t271 + t268 * s13 + (-t255 + t262) * s26 + (-s12 + s1
     &6 - s34 + t264 + s45) * s45 + (t222 + t266 + s46) * s46 + (-t18 * 
     &t70 + t267 + t79) * s56 - t269 + t270
      t273 = t16 * s16
      t274 = s16 - s34 - s13
      t275 = t16 * t274 + t262 + t43
      t276 = s25 + s26
      t277 = (s34)**2
      t278 = t259 * s25
      t279 = (-t275 - s45) * s45
      t280 = t18 * s13
      t281 = t262 - t280
      t282 = t255 - s25
      t283 = -t10 * t70 + t16 * t222 - s12 + t43
      t268 = t268 * s25
      t284 = t282 * s34
      t285 = (t16 * (s25 + s16 - s34 - s13 + s45) + s12 + t43 + s46) * s
     &46
      t286 = t255 * s16
      t287 = t18 * t271
      t288 = s26 * t281 + s56 * t283 + t260 - t268 - t279 + t284 + t285 
     &- t286 + t287
      t289 = s16 + s26 + s56
      t290 = t18 * s34
      t291 = t289 - t290
      t292 = t16 * s34
      t293 = t289 - t292
      t229 = t229 - t249
      t294 = s13 * t291 + t229 * t293 + t54 * t293
      t295 = intHs16s25s26s34s56x1112D2eps1()
      t296 = t16 * (s13 - s45 - s46 - s56)
      t297 = -t296 + s25
      t298 = intHL0s25s26s34s56x1220D4eps1()
      t299 = t48 + t292
      t300 = s25 + s26 + s34 + s56
      t301 = t48 + t249
      t302 = s13 * t299 - t300 * t301
      t303 = intHL0s25s260s56x1022D4eps1()
      t304 = s26 + s56
      t305 = t60 + t304
      t306 = s45 + s46 + s56
      t307 = t16 * t306
      t308 = s13 * t305
      t309 = t304 * t306
      t310 = (-t307 - s26 - s25) * s25 + t308 - t309
      t311 = intHs160s26s34s56x1013D4eps1()
      t312 = intHL0s25s26s34s56x1210D2eps1()
      t313 = intHL0s25s26s34s56x1310D4eps1()
      t314 = intHs160s26s34s56x1012D2eps1()
      t315 = intHs16s25s26s34s56x1111D2eps0()
      t316 = intHL0s25s260s56x1012D2eps0()
      t317 = intHL0s25s260s56x1013D4eps0()
      t318 = intHL0s25s260s56x1021D2eps0()
      t319 = intHs160000x0111D0eps0()
      t320 = intHL0s25s26s34s56x1130D4eps0()
      t321 = s13 - s45 - s56
      t322 = -t16 * t321 + s25 + s46
      t323 = intHLs160000x0211D2eps0()
      t324 = -s16 + s25 + s34
      t325 = t16 * (s12 + s13 - s56)
      t326 = t18 * s46
      t327 = t10 * s45
      t328 = -t325 + t324 + t327 + t326
      t329 = t10 * t248
      t325 = -t325 + t324 + t329
      t330 = intHs160000x0121D2eps0()
      t331 = s25 - s26 - s16 + s34
      t332 = s13 - s56
      t333 = t18 * t332
      t334 = -t333 - t54 + t331 + t329
      t335 = intHLs16s25s26s34s56x1211D4eps0()
      t336 = s25 - s26
      t337 = s16 - s34
      t338 = t16 * t337
      t339 = t42 * s46
      t340 = t41 * s45
      t341 = t339 - t333 - t338 + t336 + t340 - t61
      t262 = t16 * t262
      t342 = t41 * t248
      t343 = -t333 + t262 + s26 - s16 + s34 + t342
      t344 = intHLs16s25s26s34s56x1121D4eps0()
      t345 = s46 + s56
      t346 = t10 * s13
      t340 = t345 * t42 + s26 + t21 + t340 - t346
      t347 = s16 + s26 - s34
      t348 = s56 * t42 + t262 + t342 - t346 - t347
      t339 = -t321 * t41 + t21 + t339
      t349 = s13 + s16
      t350 = t54 + t349
      t351 = -s16 + s25
      t352 = -t351 - t292
      t353 = -t351 - t290
      t354 = t16 * t352 * t248 - (t54 + s13 + s16) * s16 + s25 * t350 + 
     &s26 * t351 - (t16 * (s25 - s12 - s16) - s26 - t280 + s34) * s34 + 
     &s56 * t353
      t355 = (t350 + s26 - s34) * s25
      t356 = (s46)**2
      t357 = s25 * s45
      t358 = t16 * (-t357 + t356) + (-t258 + s26) * s46 - (-t326 + s25) 
     &* s56 + t355
      t359 = intHLs16s25s26s34s56x1131D4eps1()
      t360 = (s26)**2
      t361 = s12 * t48
      t362 = s26 * s34
      t307 = t307 + s25
      t363 = s13 * (t106 + s25)
      t364 = t48 * t307
      t365 = -t364 + t363
      t366 = t16 * t33
      t367 = s25 + t366
      t368 = t54 * s13
      t369 = -t60 + t33
      t370 = s13 * (-t21 + t366) - t307 * t369
      t371 = s25 + s56
      t372 = t16 * t304
      t373 = t307 * (t372 + t21) - t255 * (t16 * t371 + s26)
      t374 = (t10 * t33 + t54) * s13 - t307 * (-s25 + t366)
      t375 = intHLs160000x0121D2eps1()
      t376 = propW16 * s16
      t377 = t376 * t223
      t378 = t16 * (t377 - s25)
      t379 = t377 * t41
      t380 = t19 * (t29 * (-t18 * t345 - s26 + t280 - t327 - t60) + t379
     &) + t22 * (t18 * (s13 - s46 - s56) - s26 - t327 + t378)
      t353 = s13 * t353 - s16 * t301 + s25 * t48 + t10 * s34 * t248 + t1
     &6 * (s34 + s45 + s46) * s25 + t362 + t290 * s56
      t381 = intHLs16s25s26s34s56x1211D4eps1()
      t382 = ddreal(7)
      t383 = t18 * t337
      t384 = ddreal(8) * s45
      t385 = t10 * t332
      t386 = t382 * s46
      t256 = t16 * t256
      t387 = ddreal(8) * t248
      t388 = intHLs16s25s26s34s56x1121D4eps1()
      t389 = t42 * s13
      t390 = t10 * s25
      t391 = t382 * s56
      t392 = (s25)**2
      t393 = t389 * s25
      t342 = -t16 * (s13 * t304 - s26 * s45 + (-s26 - s45 - s56) * s56) 
     &+ t18 * (s46 * s56 + t392) + s25 * (t342 + t391 + s26) + s26 * s46
     & - t393
      t394 = t16 * t259
      t395 = t18 * s16
      t396 = s12 + s16
      t397 = t16 * t396
      t398 = t60 - t280 + s26
      t399 = s12 * s13
      t400 = t248 * t324
      t401 = t398 * s34
      t402 = (t79 - t346 + t21) * s34
      t403 = (t394 + t395 - s25) * s25
      t404 = (s12 - s13 + s16 + s25 + s26) * s26
      t405 = s13 * s16
      t406 = t16 * ((t338 - s25 + s12 + s26) * s56 - t399 + t404)
      t407 = t10 * (t400 + t405)
      t408 = s25 * s26
      t409 = s26 - s13 + s45 + s56
      t410 = s12 - s45
      t411 = s26 - s34
      t412 = -t16 * (s16 - s13 + s45) - t411
      t413 = s12 + s16 + s26 - s34
      t414 = t16 * t413
      t415 = (t414 + t326) * s56
      t416 = (-s13 + s25 + s26) * s34
      t417 = t396 * s13
      t418 = t16 * (s25 * t396 - s45 * t324 + t356 + t404 - t416 - t417)
      t419 = t346 - t21
      t420 = s26 + s16 - s34 - s13 + s45 + s46
      t254 = t254 * s25
      t331 = t248 * t331
      t274 = (t16 * (s12 + s45 + s46) + t18 * t274 + s25 + t37) * s56
      t421 = (t18 * t263 + s16 + s26 + t54) * s26
      t422 = t16 * (t331 + t399 - t254)
      t423 = t18 * (-t405 + t271)
      t424 = intHLs16s25s26s34s56x1221D4eps1()
      t425 = -s13 + s16
      t426 = s25 + s26 + s34
      t427 = s25 + s45 + s46
      t428 = t10 * t427
      t429 = t42 * s34
      t430 = t43 + t60 - s13
      t431 = t16 * t410
      t432 = -t431 + s25 - s16 + s34 - s13 + s46
      t433 = t426 * s46
      t434 = (t394 + s16) * s25
      t435 = (-t54 - s25 - s16 - s13 - s26) * s26
      t436 = t16 * (s45 * t426 - t399)
      t437 = (-t54 - t425 - s26) * s26
      t350 = (-t350 + s25) * s25
      t438 = t60 + s26
      t439 = (-s12 - s16 - s25 - s26) * s26
      t258 = t16 * (s34 * t276 + s45 * t438 + (s25 - s12 - s13 + s45 - s
     &46 + s56) * s56 - t356 + t439) - (t61 + t280 + t273) * s25 + (t258
     & + s25) * s46
      t440 = t16 * t420 + s25
      t441 = s26 + s16
      t264 = (-t264 - s25) * s34 + t16 * ((t16 * t441 + s12 - s13 + s25 
     &- s34 + s45 + s46 + s56) * s56 + t248 * (t411 + t273) + t270 - t39
     &9 + t404) - t346 * s16
      t270 = t425 * s25
      t404 = -t282 * s26
      t442 = t248 * t347
      t443 = t442 - t417 + t271
      t444 = t60 - t389 + s16
      t445 = s13 * t382 - s26 - t390
      t446 = s26 + s45 + s46
      t447 = t18 * t425
      t448 = t10 * t446
      t449 = t382 * s34
      t450 = t10 * t259
      t451 = t16 * (s26 + s16 - s13 + s45 + s46)
      t452 = t248 * (t441 - t292)
      t453 = t16 * (t452 - t405 + t271 - t277) + (t349 + s26) * s25 + (t
     &43 + t450 + t273) * s34 + (t451 - t290 + s25) * s56
      t425 = t60 * t425
      t454 = t16 * t282
      t455 = s25 - s34
      t394 = t16 * (t442 - t405 + t271) + s25 * t349 + (t441 + t394 - s3
     &4) * s34 + (t451 + t455) * s56 + t404
      t442 = s26 + s34
      t451 = s25 * t248
      t456 = t10 * t451 + t16 * (s34 * t306 + t392) + s13 * (t338 - t21)
     & - s16 * t307 + t442 * s25 + t21 * s56
      t457 = intHLs16s25s26s34s56x1141D6eps0()
      t458 = t427 * s16
      t269 = t16 * (t405 + t269) - (s26 + t397 + t280 - s25) * s25 + (-t
     &338 + t390) * s45 - (t338 - t390 + s26) * s46 + (-t338 + t21 + s46
     &) * s56
      t459 = t10 * t417 - t287 - t425 + (-t444 - s26) * s26 + (-t389 + t
     &21 + s26) * s34 + (-t60 - t447 + t429 - t448) * s56 + t248 * (s34 
     &* t41 - t10 * t441)
      t336 = t16 * (t399 + t392) - t287 - (t395 + s13) * s25 - (t82 - t3
     &89 + s26) * s26 - s34 * t445 + (-s16 * t42 + s25 + t280 - t448 + t
     &449) * s56 + t248 * (-s16 * t41 + t10 * t336 + ddreal(8) * s34) + 
     &t389 * s16
      t460 = t80 * s34
      t461 = -t16 * ((s25 + s26 - s16) * s45 - t254 + t271 + t405 + t439
     &) + (-t276 + t273) * s46 + (-t16 * (s25 - s12 - s16 - s13 + s45) -
     & s46) * s56 - t460
      t462 = intHLs16s25s26s34s56x1321D6eps0()
      t463 = s25 + s16 - s13
      t464 = t10 * s34
      t262 = -t16 * (t248 * (t276 + t292) + t399) + (t54 + s16) * s25 + 
     &(t54 + t463 + s26) * s26 + (-t262 + t280 + s16 - s34) * s34 + (-s2
     &5 + s16 + s13 + t225 - t464 - s56) * s56 - t405
      t292 = -t16 * (t248 * (t441 - t290) - t417) + (-t82 + t280 - s26) 
     &* s26 + (-t16 * t446 - s56 - t463 + t464) * s56 - t270 + t402
      t441 = intHLs16s25s26s34s56x1131D4eps0()
      t446 = t16 * t276
      t463 = intHLs16s25s26s34s56x1221D4eps0()
      t465 = intHLs16s25s26s34s56x1231D6eps0()
      t466 = t16 * (t248 * (t446 + s34) + t399) + t287 + (t16 * t263 + s
     &26) * s34 + (t18 * (s34 - s13) - s16 - t256 + t428) * s56 + t350 +
     & t405 + t437
      t467 = t382 * s25
      t263 = t10 * (t60 + t442) * t248 + t16 * (s16 * t70 + t392 + t399)
     & + t287 - (t450 + t395) * s25 + (t10 * t263 + s26) * s34 + (-t18 *
     & t349 - t256 + t329 + t429 + t467) * s56 + t435
      t349 = intHs160000x0211D2eps0()
      t442 = intHL0s25s260s56x1031D4eps0()
      t468 = s26 - s56
      t469 = intHL0s25s260s56x1022D4eps0()
      t470 = -t16 * (s12 + s26 + s13 - s45 - s46 - s56) - t337
      t471 = intHL0s25s260s56x1021D2eps1()
      t472 = t25 * s26
      t473 = epinv * s56
      t474 = intHs160000x0111D0eps1()
      t475 = intHL0s25s26s34s56x1130D4eps1()
      t476 = intHs16s25s26s34s56x1113D4eps0()
      t477 = intHLs160000x0211D2eps1()
      t478 = intHs160000x0121D2eps1()
      t479 = intHLs160000x0121D2eps0()
      t378 = t10 * t321 - t326 + t378
      t326 = t29 * (-t10 * (s45 + s56) - t326 - t60 + t346) + t379
      t480 = t18 * epinv
      t481 = intHLs16s25s26s34s56x1231D6eps1()
      t482 = t25 * s34
      t483 = -t482 + s26 + s16 + s56
      t484 = t103 * s34
      t485 = t484 + t33
      t486 = intHs160000x0211D2eps1()
      t487 = t480 - t10
      t488 = intHL0s25s260s56x1031D4eps1()
      t489 = intHs16s25s26s34s56x1112D2eps0()
      t490 = intHLs16s25s26s34s56x1311D4eps1()
      t491 = intHL0s25s26s34s56x1220D4eps0()
      t492 = intHL0s25s26s34s56x1120D2eps1()
      t493 = epinv * s34
      t494 = t493 - s25 - s26 - s56
      t495 = intHLs16s25s26s34s56x1321D6eps1()
      t496 = intHs16s25s26s34s56x1211D2eps0()
      t497 = epinv * t250
      t498 = t79 + t222
      t499 = propW16 * t223
      t500 = t499 * t48 * t33
      t501 = ddreal(1) / t48
      t502 = ddreal(1) / (ecossin)**2
      t503 = t298 * t25 * t302 * t501
      t247 = propW34 * (t135 * (-t101 * (t10 * (t452 - t417) + t287 + ep
     &inv * t453 + (t444 + s26) * s26 + s34 * t445 + (t60 + t447 - t449 
     &+ t448) * s56 + t425) - t118 * t453 - t119 * (t16 * (s12 * t247 + 
     &s26 * t304 - s34 * t427 - t451 + t458) + t18 * (s56 * t337 + t408)
     & + s13 * (-t383 - t79 + s25) - s25 * s56 + s26 * t337) + t128 * (e
     &pinv * t264 + s56 * t440 + t16 * t443 + t270 + t284 + t404) - t13 
     &* (-t398 * s34 + t274 + t421 - t422 + t423)) - t250 * t496) * t24
      t270 = t208 * propW34
      t258 = t270 * (-t104 * (t16 * ((t265 + s46) * s46 - t357) - (-s46 
     &* t16 + s25) * s56 + t355) - t105 * (t10 * (t248 * t438 + t271 + t
     &392) - t393 + epinv * t258 + (-t255 + t21) * s26 + (ddreal(9) * s2
     &5 + t10 * (s26 - s13 + s45 + s46)) * s56) - t109 * t461 - t120 * t
     &263 + t126 * (-s46 * t412 + t415 + t418) - t155 * t258 - t424 * (t
     &10 * t426 * t248 + t287 + t368 + epinv * (t436 + s34 * t430 + (t43
     &2 + s56) * s56 + t433 - t434 + t435) + (-t346 + t21 + s26) * s34 +
     & (-t256 + t429 + t428 - t280 - s16) * s56 + t405 + t437 + t350) + 
     &t47 * (epinv * (-t16 * (s13 * (s25 - s16 + s56) + (s25 + s16 - s34
     & - s45) * s26 + (-s25 + s16 - s45) * s56 - t271 - t357 + t360 + t3
     &61 + t458) + s25 * s34 + s46 * t48) - t363 + t364) + t463 * (-t436
     & - s34 * t430 + (-t432 - s56) * s56 - t433 + t434 - t435) + t503)
      t256 = -t121 * (-t406 + epinv * (t16 * (t400 + t399 - t360) - (t39
     &7 + s13) * s25 - (t16 * (s12 - s13) + s16 + t21) * s26 + (-t383 - 
     &t83 + s25) * s56 + t401 + t280 * s16) + t402 - t403 + t407) - t131
     & * (epinv * t348 + s16 + s26 - s34 - t21 - t387 + t389 - t391 + t5
     &4) - t133 * (epinv * t343 + s16 - s34 - t21 + t256 + t385 - t387) 
     &- t143 * (t10 * (t331 + t405 - t271) + t368 + epinv * (-t423 + t42
     &2 - t421 + t401 - t274) - (t54 + t395 + s13 - s25) * s25 + s26 * t
     &419 - t419 * s34 + (-t10 * t420 + s25) * s56) - t250 * t319 - t348
     & * t7
      t265 = epinv * t375
      t274 = t502 * propW34
      t45 = t274 * (t135 * (t24 * (-t3 * t325 - t315 * (s13 + s26 - s56)
     & - t330 * t334) + t27 * t256 + t31 * (-t113 * t322 - t116 * (epinv
     & * t339 + ddreal(8) * t321 - t386 - t390) - t123 * (-t406 + epinv 
     &* (s46 * t412 - t415 - t418) + t402 - t403 + t407) + t2 * (t54 + s
     &25 + s26 + s16 - s34 - s13 - s46 + s56) - t28 * t339 + t49 * (-epi
     &nv * (t16 * (s25 * t410 + s46 * t409 + t356) + t224 * s25) + s25 *
     & t250) + t323 * t328 - t335 * t341 - t344 * t340 - t381 * (epinv *
     & t341 - s25 + t383 - t384 + t385 - t386 + t45 + t79) - t388 * (epi
     &nv * t340 - t345 * t382 - s26 - t384 + t389 - t390)) - t265 * t380
     &) - t31 * (t316 + t318) * t501 * t64)
      t224 = t240 + t241
      t256 = t244 + t245
      t284 = t501 * t64
      t321 = (t242 - t243) * t501
      t331 = t24 * s34
      t339 = t274 * t135
      t257 = t339 * (t158 * (t27 * (t261 * (t182 - t110 - t111) + t171 *
     & t294) + (-t231 * t238 - t232 * t95) * t24 * t17 - t24 * (t251 + t
     &253) * t252 * s34) + t27 * (t163 * t336 - t343 * t6) + t31 * ((t39
     & * (epinv * t322 + s25 - t296) - t284 * t256) * s34 + t501 * (-t30
     &2 * t491 - t492 * (s13 * (-t484 + t48) + t301 * t494) + t17 * t64 
     &* t224) + t321 * t64 * t25) + t27 * (-t16 * (s16 * t259 - (s12 + s
     &13 + s16 + s26) * s34 - t271) - (t18 * t259 + s26 + t273) * s26 - 
     &(-t267 - t266 - s46) * s46 + (-t10 * t257 - t18 * (s12 - s46) + t2
     &76) * s56 + t260 - (s16)**2 - t277 - t278 - t279) * (-t161 + t68) 
     &* (t158)**2 - t331 * t297 * (t129 * t295 + t489))
      t259 = t25 * intHL0s25s260s56x1012D2eps1()
      t267 = t132 * t487
      t127 = t274 * (t135 * (t27 * (-t267 * t297 - t474 * t497) - epinv 
     &* (t127 * t325 + t334 * t478) * t24) + t31 * (t135 * (t167 * t263 
     &+ t174 * t322) + t501 * (t135 * (-t471 * (-t16 * (s13 * (-t99 - s2
     &6 - s56) + t309) - s25 * (t99 + t16 * (s45 * t25 + s46 * t25) + t7
     &1 - t473 + t472)) + t64 * (t312 * t482 - t234 + t235 + t236 + t58)
     &) + t259 * t64) + t135 * (-t114 * t322 + t328 * t477) * epinv) + t
     &135 * (t24 * (t252 * t314 * t482 - t228 * t231 - t232 * t67) + t27
     & * (-t294 * t5 - t4 * (s13 * t485 + t229 * t483 + t483 * t54) - t9
     &8 * (-t287 + t286 + epinv * t261 - s26 * t281 - t282 * s34 + (-t27
     &5 - s45) * s45 - s56 * t283 - t260 + t268 - t285))) * t158)
      t45 = t127 + t257 + t45 + t502 * (-t135 * (t479 * propW34 * t380 +
     & t1 * propW34 * (t19 * t326 + t22 * t378) + t115 * (t19 * propW34 
     &* (epinv * t326 - t379) + t107 * propW34 * (t378 * t480 - t379)) +
     & epinv2 * (t19 * (t16 * propW34 * (t16 * (t248 * (t43 + t60 + t337
     &) - t405) + t287 - (t450 + s16) * s25 - (s13 * t41 + s26 + t61 + t
     &82) * s26 + (t255 + t276) * s34 - (t41 * t70 - t21 + t337 + t61 - 
     &t79) * s56 + t11 * (s56 * t413 + t254 - t400 - t405 - t416 - t439)
     &) + t500 * t41 * propW34) + t136 * propW34 * s25 * (-s13 * t337 - 
     &(t255 + s12) * s26 + (s25 - s12 + s26 - t72 + s56) * s56 + t248 * 
     &t498 - t278 + t500)) * t501 * t158) + t247 + (-t135 * (t108 * t336
     & + t140 * t264) + t138 * (t497 - s26 - s16 + s34 - s13 + s56 - t22
     &5)) * t27 * propW34 + t258)
      t70 = t135 * t24
      t82 = t70 * (-t141 * (s56 * t16 * t440 + epinv * t394 - s26 * t454
     & + s34 * t454 + t10 * t443 + t425) - t297 * t349 - t394 * t9)
      t127 = t24 * t158
      t136 = epinv * t486
      t225 = t239 * t158
      t229 = t274 * (t135 * ((-t225 * t24 * t231 - t476 * t27 * t470 - t
     &490 * t31 * (epinv * t233 - s16 + s25 - s26 + s34 + t329 - t333 - 
     &t54)) * s34 + t24 * (-t158 * (t261 * t32 + t272 * t96 + t288 * t97
     &) - t136 * t297) + t25 * (t27 * (t292 * t94 + t353 * t87 + t456 * 
     &t85 + t459 * t93) + t31 * (-t250 * t475 - t262 * t495 + t269 * t92
     & - t297 * t488 + t342 * t90 - t354 * t81 + t466 * t481) + t127 * (
     &t261 * t66 + t272 * t52 + t288 * t53)) + t320 * t31 * t250) + t284
     & * t31 * (s34 * t135 * t246 + t317))
      t232 = t78 * t25
      t247 = s34 * intHLs16s25s26s34s56x1311D4eps0()
      t252 = t25 * intHL0s25s260s56x1013D4eps1()
      t254 = t311 * t158
      t64 = t274 * (t31 * (t135 * (t501 * (-t313 * t482 * t64 - t310 * t
     &469) + t232 * t373 - t247 * t233) - t252 * t284) + t135 * (t25 * (
     &t139 * t226 + t370 * t74 + t374 * t76) - t211 * (epinv * t470 + s2
     &6 + t329 - t333 + t60)) * t27 + t254 * t331 * t209 * t231)
      t64 = t64 + t229 + t274 * (t27 * (-t135 * (t134 * t226 + t146 * t3
     &70 + t147 * t374 + t150 * t456 + t153 * t353 + t292 * t63 + t459 *
     & t59) - t14 * t297 - t15 * (epinv * t297 - s16 + s25 - s26 + s34 +
     & t329 - t333 - t54)) + t208 * (-t148 * t373 + t149 * t354 + t154 *
     & t358 - t157 * t342 + t297 * t442 - t359 * (-epinv * (t16 * (-s45 
     &* t48 + t361 - t362) + (s16 - s34 - s46) * s25 + (s25 + s16 - s46)
     & * s26 + (-s25 + s16 - s46) * s56 - t271 + t308 + t360) - t16 * t3
     &65) - t57 * t269 - t69 * ((t60 + t327 - t280) * s25 - epinv * t358
     & + (t390 - s26) * s46 + (t21 + s46) * s56 + t408) + t462 * t262 + 
     &t441 * (-s46 * t276 - (-t431 + s25 - s16 - s13 + s46 + s56) * s56 
     &+ t434 - t435 - t460 - t446 * s45) - t465 * t466 + t25 * t303 * t3
     &10 * t501) + t82)
      t82 = t25 * intHLs16s25s26s34s56x1411D6eps1()
      t139 = t25 * intHLs16s25s26s34s56x1141D6eps1()
      t229 = t139 - t457
      t156 = -t183 + t156
      t183 = -t187 + t152
      t145 = -t185 + t145
      t185 = (t82 - t227) * s34
      t231 = t25 * intHLs160000x0122D4eps1()
      t233 = intHLs160000x0111D2eps1() * t17
      t255 = -t25 * intHLs160000x0113D4eps1() + intHLs160000x0113D4eps0(
     &)
      t257 = s25 * t129
      t258 = epinv * t33
      t260 = -s16 + s25 - s26 + s34 - s56
      t261 = epinv * s25
      t195 = t195 * t260
      t262 = t120 * t300
      t263 = s34 * t489
      t264 = (t195 - t228) * t158
      t268 = t305 * t469
      t269 = t471 * (t99 + s26 + s56)
      t271 = t116 * t487
      t272 = t96 * t260
      t275 = t98 * (t261 + s26 + s16 - s34 + s56)
      t251 = (t251 + t253) * s34
      t253 = t25 * t131
      t277 = epinv * t478
      t222 = t135 * (t24 * (t253 - t267) + t27 * (-t140 * t351 - t277)) 
     &+ t31 * (t135 * (t105 * (t257 + s26 + s56) - t109 * t351 - t25 * t
     &388 - t501 * (t268 + t269) + t344 - t265 + t271) + t501 * (t316 + 
     &t318) + t155) + t135 * (t19 * (-t12 * t330 + (t315 + t7) * t12 - t
     &29 * t479) + t27 * (-t108 * t293 - t263)) + t31 * (t135 * (t123 * 
     &t324 - t33 * t463 - t424 * (t258 - s34) + t501 * (-t58 + t234 - t2
     &35 - t236) + t262) - t259 * t501) + t70 * (-t101 * (-epinv * t352 
     &+ s34) + t118 * t352 - t119 * t324 + t13 * t369 + t143 * (epinv * 
     &t369 + s16 + s26 - s34 + s56) + t264) + t135 * (t24 * (-s34 * t5 -
     & t17 * t238 - t272) + t27 * (-t251 - t275) - epinv2 * (t19 * ((pro
     &pZ25 * t324 - t10) * s25 - t304 * t41 - t338) - t22 * (t372 + t222
     &)) * t501) * t158
      t267 = t303 * t305
      t278 = t17 * t501
      t279 = t100 * t158
      t280 = t295 * s34
      t281 = t121 * t129
      t282 = t24 * t67
      t283 = t182 * t27
      t284 = t25 * t501
      t114 = t114 * epinv
      t285 = t220 * t135
      t286 = (t135 * (t27 * ((-epinv * t4 + t25 * t314) * s34 - t68) - t
     &282) + t283) * t158
      t287 = epinv * t359
      t225 = t19 * (t12 * t14 - t135 * t168) + t25 * (t24 * (t135 * ((t8
     &7 + t254) * s34 + t260 * t74 + t324 * t85) + t158 * t66) + t208 * 
     &((t313 * t501 - t81) * s34 - t305 * t90 - t324 * t92 - t48 * t78) 
     &+ t76 * t27 * t33 * t135) + t31 * (t441 + t69) + t135 * (t19 * (-t
     &12 * (t349 + t97) - t29 * t442) + t24 * ((-t153 - t476) * s34 - t1
     &50 * t324 + t260 * (-t9 - t146)) + t27 * (-t147 * t33 - t136 + t19
     &3) + t31 * ((t490 + t149) * s34 + t148 * t48 + t157 * t305 + t25 *
     & t488 + t324 * t57)) + t31 * (t252 - t317) * t501 + t287 * t31 - t
     &15 * t129 * t24 + t135 * (-t24 * (t129 * t151 + t225) - t31 * t246
     & * t501) * s34 - t135 * (t141 * t24 * t260 + t115 * t31) * epinv -
     & t127 * t32
      t254 = -t188 + t56
      t260 = t254 * t24
      t288 = t28 * t31
      t292 = ddreal(9) * t73
      t294 = t376 * t41
      t295 = t11 + t294 + t16
      t301 = t376 * t16
      t302 = ddreal(1) + t301
      t303 = t41 * propW16
      t304 = s25 + s26 - s34 + s56
      t308 = t1 * t295
      t309 = epinv2 * (t19 * (t303 * t33 - t10 + t11) + t22 * (propW16 *
     & t366 - ddreal(1))) * t158
      t310 = t171 - t5
      t95 = (-t95 - t238) * t17
      t238 = t115 * ((t19 * t295 + t22 * t302) * epinv - t294 * t125)
      t295 = t1 * t302
      t302 = ((-t265 - t231 - t479 + t237) * s16 + t233 + t91) * t135 * 
     &t125 * propW16 + t73 * (-t191 + t97) - t192
      t239 = -t25 * t311 + t239
      t311 = t25 * t81 - t149
      t322 = t25 * t90
      t325 = t322 - t157
      t326 = t25 * t92
      t327 = t326 - t57
      t328 = -t221 + t32
      t193 = t135 * (t31 * (t311 * s34 + t154 * t455 + t305 * t325 + t32
     &4 * t327 + t48 * (t232 - t148) - t69 * (t493 + t99) - t28) - t193 
     &* t27 + t331 * (t239 * t158 + t151)) + t127 * t328 + (t141 + t143 
     &+ t15) * t24 * epinv + t26 * (t135 * t97 - t13 - t14 - t9)
      t305 = t376 * t125
      t329 = t305 * t135 * t255
      t136 = t16 * t193 + t18 * t73 * (t295 + t2 + t7 + t8 - t330 - t349
     &) - t41 * t302 - ddreal(12) * t329 - t135 * (t24 * (-t281 * t48 + 
     &t253) + t31 * (-t167 * t300 + t113) - t238) - t158 * (t27 * (t135 
     &* (-t251 + t161 - t275) - t110 - t111) + t70 * (s34 * t310 - t272 
     &+ t95)) - t135 * (t19 * (t12 * (-t330 - t349) - (-t7 - t8) * t12 -
     & t308) + t24 * (-t101 * (epinv * t48 + s34) + t132 * t25 + (t195 -
     & t228) * t158 + t48 * (-t118 - t119)) + t27 * (-t140 * t304 - t136
     & - t280) + t31 * (-s34 * t104 + t105 * (t493 + t99 + s26 + s56) + 
     &t123 * (epinv * t289 + s25 + s26 + s56) + t126 * t289 - t155 * t45
     &5 + t262) - t309) - t135 * (t27 * (-t128 * t304 - t478) + t31 * (-
     &(t49 + t39) * s34 - t166)) * epinv - t286 + t285 - t208 * (-t176 +
     & t177 - t174 + t114)
      t166 = t16 + epinv
      t193 = t16 * (s12 + s26 + s16) - t290
      t195 = -t414 + s25
      t228 = s12 + s26 + s16 - s34 + s56
      t253 = t10 * t228
      t228 = t16 * t228
      t262 = -t228 + s25
      t289 = t43 + t397 + t371 - t464
      t290 = -s13 + s25 + s26 + s45 + s56
      t302 = t10 * t396
      t304 = t411 * t42 + t302 + t60 + t71
      t333 = t253 + s25
      t334 = t376 * t10
      t336 = ddreal(1) + t334
      t338 = t42 * s26
      t340 = t382 * t411 + t396 * t41 + s25 + t71
      t341 = t16 * t409 + t21
      t342 = t10 * s56
      t266 = t266 + t342 + t21
      t343 = -t16 * (s16 - s34 + s13 - s45 - s56) + s25 - t61
      t345 = ddreal(12) * t376
      t348 = t11 + t345 + t16
      t350 = t15 * t166
      t205 = (-t205 - t206) * t25
      t176 = t208 * (-t176 - t174 + t84 + t210 + t114)
      t95 = t158 * (t135 * (t27 * (t314 * t482 - t68) - t282) + t283) + 
     &t135 * (t24 * (t158 * (t95 - t272) + t25 * (t195 * t85 + t262 * t7
     &4)) - t238) + t27 * (-t158 * (t135 * t275 + t110 + t111) - t251 * 
     &t158 * t135 + t209 * (t193 * t87 + t289 * t94 + t304 * t93 + t333 
     &* t76 + t279)) + t135 * (t19 * (-t12 * (t349 + t330 + t134) - t308
     &) + t31 * (t441 * t48 + t69 * (t103 * t290 + s25) + t247) - t309) 
     &+ t24 * (t135 * (-t140 * t396 - t141 * (t102 - t253) - t143 * (t26
     &1 - t228) - t146 * t262 - t150 * t195 + t264) - t350) + t135 * (-(
     &t478 + t486) * epinv - t147 * t333 - t153 * t193 - t289 * t63 - t3
     &04 * t59 + t217 - t280) * t27 + t208 * (epinv * ((-t39 + t490) * s
     &34 + t48 * (-t47 + t359)) + t113 + t205) + t26 * t9 + t176 + t26 *
     & (t13 + t14)
      t176 = t300 * t463
      t193 = epinv * t477
      t195 = t300 * t424
      t29 = t135 * (t27 * (t158 * (-t291 * t310 + t4 * t485) - t163 * t3
     &40) + t31 * ((t341 * t49 - t195) * epinv + t38 * (-t120 + t167)) +
     & t220 + t180 * t24) + t137 * t6 + t27 * (t135 * (t108 * t340 + t17
     &0) + t138) + t208 * (t104 * t341 - t105 * (epinv * t266 + s25 + t7
     &9) - t123 * (epinv * t343 + s25 + t79) - t126 * t343 - t155 * t266
     & + t25 * (t381 + t388) - t271 - t335 - t344 - t176 + t193) + t265 
     &* (t19 * t348 + t22 * t336) * t135 + t70 * (t101 * (t16 * (epinv *
     & s26 - s25) + t261 - t302 + t449 - t71 - t338) + t121 * (-t10 * (s
     &12 + s16 - s34) - t16 * t472 + t257) + t128 * (t16 * (t129 * t396 
     &+ s26 - s34 + s56) + s25) + t160 * t487 + t162) + (t12 * t3 + t29 
     &* t323 + t348 * t479) * t135 * t19
      t162 = t26 * t173
      t206 = (-t231 + t237) * s16 + t233 + t91
      t220 = (-t123 - t105 + t49) * epinv
      t238 = epinv * t69 + t154
      t251 = t43 + t54 + t390 + s56
      t253 = s25 - s26 + s56
      t257 = s16 + s56
      t262 = t106 + s25 - s13 + s45 + s46
      t264 = t71 + t276 - t72
      t266 = t129 * s12
      t271 = t16 * (t266 + s25 + s26)
      t76 = t25 * t76
      t272 = t25 * t481
      t275 = t490 * s34
      t59 = (t25 * t93 - t59) * t62
      t62 = (t25 * t94 - t63) * t65
      t15 = t19 * (t12 * t9 + t135 * (-t12 * t7 - t219)) + t24 * (t135 *
     & (-t170 - t218 - t216 + t217) - t158 * t328 - t202 * t158 * t135 -
     & t203 * t158 * t135) + t27 * (t135 * (t59 - t211 + t62) - t146 - t
     &15 - t150) + t31 * t441 + t27 * (t25 * (t85 + t74) + (t25 * t87 - 
     &t153) * t135 * t324) + t208 * (t154 * t264 + t155 * t262 - t25 * t
     &475 + t311 * t324 + t359 * (t102 - t37) + t69 * (epinv * t264 - s2
     &5 + s26 - s56) + (t25 * (t90 + t92) - t157 - t57) * t253 + (t272 -
     & t465) * (t54 + t257 + t338 + t21) + (t25 * t495 - t462) * (t43 + 
     &s16 + s34 + s56 + t36) + t320 - t275) - t30 * t335 * t135
      t37 = epinv * t474
      t63 = t208 * ((-t165 + t477) * epinv + t463 * t55 + t113 + t28 + t
     &323 + t491 + t492)
      t12 = t137 * (-t3 - t134) + t24 * (-t129 * t138 + t135 * (-t172 + 
     &t171 - t169) - t496) + t135 * (t27 * (-t117 * t131 - t119 * t223 -
     & t121 * (epinv * t498 + t271) - t13 * (t54 + t33) - t133 * t487 - 
     &t143 * (t130 + t258 + t21) + (t163 - t108) * (t54 + s26 - s16 + s3
     &4 - s56 + t21) - t37) + t31 * (t105 * (t16 * (epinv * t262 - s26) 
     &- t21) - t109 * t253 - t117 * t381 - t120 * t251 - t123 * (epinv *
     & (-t71 + t411 + t72) + t271) - t126 * (-t71 + t411 + t230) - t17 *
     & t388 + t25 * (t251 * t34 - t116 - t298) + t424 * (-t25 * t54 - s1
     &6 - s56 - t21 + t261 - t338) + t47 * (s26 * t166 - t473 + t99) - t
     &49 * (epinv * t51 - s25)) - t277 * t24 - t30 * t344) + t135 * t19 
     &* (-t12 * t319 + (t315 + t330) * t12 - t168) + t63 + (t100 * t209 
     &- t110 - t111) * t158 * t27
      t12 = -t10 * t135 * (t31 * t276 * (t232 - t148) + (t76 - t147) * t
     &35 * t27) - t16 * t15 + t18 * (t207 * t26 + t73 * (t344 + t315 - t
     &319 + t330 + t1 + t4 + t5)) - t41 * (t135 * (t27 * (-t190 - t189) 
     &+ t31 * ((t82 - t227) * s34 - t184 * t468 + t80 * (t139 - t457))) 
     &+ t73 * (-s25 * t9 - t2 - t335 + t7)) - t135 * (t24 * (-t118 * t80
     & - t181) + t31 * (-t179 - t175 - t174) + t23 * (-t4 - t5)) - t12 -
     & (-t135 * t68 + t182) * t158 * t27 - t292 * t6
      t15 = -s12 - s13 + s25 + s45 + s46 + s56
      t23 = s12 - s13 + s25 + s45 + s46 + s56
      t30 = s12 + s13 + s16 + s26 - s34 - s45 - s46
      t35 = -s13 + s25 + s45 + s46 + s56
      t51 = t16 * t48
      t63 = t71 + t390 - s26
      t44 = t79 + t44 + t342
      t54 = s25 * t41 + t10 * (s34 + s56) - t395 - t54
      t43 = t16 * (s25 - s12 - s16 + s34) + s56 - t43
      t65 = t48 + t273
      t72 = -t71 + t346 - t249 - s26
      t82 = t10 * t411 + t16 * t257 + s25
      t85 = t16 * (s26 + s45 + s46 + s56) + s25 - t346
      t87 = -t16 * t347 + s25
      t83 = -t83 - t395 + s25 + t464
      t71 = -t71 + t124 - t249 + t450
      t93 = ddreal(1) + t103
      t94 = t19 * (-t11 + t345 - t16) + t22 * (-ddreal(1) + t334)
      t99 = t25 * t158
      t72 = t339 * (-t24 * (t13 * t33 + t140 * (s12 - s16 - s26 + s34 - 
     &s56) + t96 * t35 * t158 + t280 * epinv) + t27 * (t25 * (t200 + t19
     &6 + t197 + t198 + t199) - t212 - t213 - t211 * t93) + t31 * (-t154
     & * t72 - t155 * (t16 * t332 - t248) + t25 * (t81 * t83 + t321) - t
     &43 * t57 - t69 * (epinv * t72 - s26 + s56 + t21) + t278 * t224 - t
     &492 * t494 * t501))
      t43 = t339 * (t25 * (t31 * (t43 * t92 + t44 * t78 + t63 * t90) + t
     &201 * t27) + t127 * (-t15 * t32 - t23 * t97) + t31 * (-s34 * t256 
     &+ t300 * t491 + t268 + t269) * t501)
      t11 = t43 + t72 + t274 * (t24 * (t135 * (-t141 * (epinv * t367 - t
     &60 - t61) - t367 * t9 - t263 + t99 * (t15 * t66 + t23 * t53 + t35 
     &* t52)) - t496) + t27 * (-t135 * (t143 * (-t258 + t36) + t86 + t89
     &) - t14) - epinv2 * (t19 * (t303 * t48 + t11 + t16) + t22 * (propW
     &16 * t51 + ddreal(1))) * t135 * t501 + t208 * (-t148 * t44 - t149 
     &* t83 - t157 * t63 - t2 - t28 + t284 * (s34 * t312 - t298 * t300 -
     & t267)) + t135 * (t265 + t479) * (t19 * (t11 - t294 + t16) + t22 *
     & (ddreal(1) - t301))) + t274 * (t27 * (-t135 * (t37 + t319 + t214 
     &+ t215) - t350) + t31 * (t135 * (t193 + t323) + t501 * (t135 * (t5
     &8 - t234 + t235 + t236) + t259 - t316 - t318) - t135 * (t381 + t38
     &8) * t117) + t70 * (-(t478 + t486) * epinv - t134 + t217 - t330 - 
     &t349))
      t19 = t339 * (t24 * (t118 * t82 + t158 * (t159 + t67)) + t27 * (t1
     &19 * t87 + t121 * (epinv * t87 + t21 + t61 + t79) + t158 * (t15 * 
     &(t110 + t111) - t98 * (t25 * t35 + t266) + t68 * t30 * t158) + t46
     & * (-t163 + t108)) + t31 * ((t17 * t39 + t122) * s34 + t105 * (t10
     & * (epinv * s13 + s56) - t16 * (epinv * (t106 + t248) - s26) + t46
     &7) + t109 * t65 + t123 * (epinv * t71 - t10 * t337 + t21 - t79) + 
     &t126 * t71 + t49 * (epinv * t85 - t60) - t167 * t54))
      t15 = t19 + t274 * (t135 * (t1 * t94 + t115 * (epinv * t94 - ddrea
     &l(12) * t305) + t24 * (t101 * (epinv * t82 - s26 - s34 + s56 - t60
     & - t61) + t128 * (-t16 * (-t266 + t258) + s25) - t171 + t4 + t5)) 
     &+ t135 * (t133 * (epinv * t10 - t42) + t487 * (t131 + t132) - t99 
     &* (t112 * t15 + t279 * t30)) * t27 + t208 * (t104 * t85 + t120 * t
     &54 + t47 * (epinv * t65 - t51) - t113 - t114 + t177) + t138 * t93 
     &* t24)
      t11 = t10 * t274 * (t27 * (s34 * t135 * t183 + t135 * t145 * t55 -
     & t188 + t56) + t208 * t156 * t253) + t16 * t339 * (t24 * (t170 + t
     &3) + t27 * (t7 + t8)) - ddreal(4) / ddreal(3) * t11 + ddreal(8) / 
     &ddreal(3) * t274 * (t31 * (t135 * (t25 * (-t299 * t495 + t475 + t4
     &88) + t299 * t462 + t195 - t320 + t335 + t344 - t442) + t501 * ((t
     &25 * t313 - t246) * t135 * s34 - t317 + t252) + t135 * (-t272 + t4
     &65) * (t51 + s34)) + t135 * (s34 * t476 + t6) * t27) - ddreal(16) 
     &* t305 * t339 * t255 + ddreal(16) / ddreal(3) * t270 * t502 * (t35
     &9 * t48 + t275) + ddreal(2) / ddreal(3) * t15 - ddreal(8) * t339 *
     & (t31 * (t229 * t48 + t185) + t125 * t206 * propW16)
      t15 = t296 + s25
      t19 = epinv * t80
      t22 = t296 + t61 + s25
      t6 = -t134 + t6
      t9 = t24 * (-t86 - t89 + t216 - t217) + t31 * (-t84 + t247 - t210)
     & + t24 * ((t204 + t9 + t13 + t14) * s25 + t141 * t142 + t143 * t14
     &4 + t25 * (t199 + t200) + t75 * (t25 * t74 - t146) + t77 * (t76 - 
     &t147) + t59 + t62) + t31 * (t48 * (t441 + t287) - t69 * (t103 * t5
     &0 + s25) - t205 + t493 * t490) + t26 * t6 + t211 * t27
      t4 = t10 * t31 * (-t154 * t50 + t276 * (t322 - t157) + t88 * (t326
     & - t57)) + t16 * t9 - t18 * t107 * t173 + t41 * (t24 * (-t254 * s2
     &5 - t189 - t190) - t107 * t6 + t192 * s25) - t24 * (t121 * (t79 - 
     &t19 + t61 + t21) + t117 * t133 + t160 * t17 - t170 + t171 - t180) 
     &- t31 * ((t477 + t375) * epinv + t105 * (-epinv * t15 + s25 + t79)
     & + t123 * (-epinv * t22 + s25 + t79) - t15 * t155 + t25 * t381 - t
     &175 - t176 + t178) + t162 - t24 * (-s25 * t138 + t101 * (t60 - t19
     & + s26 + s34 - s56 + t61) + t46 * (-t163 + t108) - t169 - t3 - t4 
     &- t5) - t31 * ((-t164 - t195) * epinv + t116 * t17 - t126 * t22 + 
     &t25 * (-t34 * t38 + t388) - t2 + t28 + t323 - t335 - t344 + t479)
      t3 = t274 * (t135 * (struc60PM * (t10 * t31 * t238 - t16 * (t31 * 
     &(t220 - t155 - t126 + t104) + (t25 * (t52 + t53 + t66) - t32 - t96
     & - t97) * t158 * t27) - t158 * (t24 * (t25 * (t98 + t112) - t110 -
     & t111) + (t161 - t68) * t158 * t27)) + struc9PM * t4) + t20 * stru
     &c10PM + (-t10 * (t135 * (-t194 * t27 + t31 * (t154 * t290 - t276 *
     & t325 - t327 * t88)) + t127 * (-s34 * t135 * t239 + t221 - t32)) -
     & t16 * t95 + t18 * (t135 * (-t162 - t288) + t73 * (t336 * t479 + t
     &2 - t3 + t323)) - ddreal(24) * t329 + t41 * (t331 * (-t187 + t151 
     &+ t152) * t135 + t135 * t186 * (t228 + s25) * t27 + t73 * (t191 + 
     &t295 - t349 - t330 - t134) + t192 + t260) + t29 - ddreal(12) * t13
     &5 * t125 * propW16 * t206 + ddreal(9) * t73 * t173) * struc5PM + t
     &12 * struc6PM)
      t1 = t274 * ((-t16 * t225 + t18 * (t135 * (-t26 * t8 - t288) + t73
     & * (t479 + t315 - t330 + t2 + t7)) + t41 * (t73 * (s25 * t14 + t1 
     &- t349 + t442 - t97) + t192 + t260) - t27 * (t135 * (t25 * (t293 *
     & t40 + t279) - t280 * t129 - epinv * t128 * t351) - t158 * (t110 +
     & t111)) - t222 - t208 * ((-t129 * t39 + t501 * (epinv * t492 + t24
     &4 + t245 - t491)) * s34 + t25 * (-t300 * t34 + t501 * (t267 - t242
     & + t243)) + t113 - epinv * t47 * t351 - t278 * (t240 + t241)) - t7
     &0 * (t158 * (s34 * t171 - t159) - t281 * t324) - t286 + t285 - t20
     &8 * (s34 * (-t122 + t284 * (-t312 + t298)) + t114) + t292 * t8) * 
     &struc2PM + (struc3PM + struc8PM) * t136 - t208 * (-t16 * t238 + t1
     &04 - t126 - t155 + t220) * struc61PM)
      result = ddreal(4) / ddreal(3) * t1 + (-t10 * t339 * t377 * t125 *
     & t255 + t16 * t274 * (t135 * (t27 * ((s13 * t367 - t307 * t33 + t3
     &68) * t145 + t183 * t250 * s34) + t31 * (t365 * t229 + t156 * (-s4
     &6 * t468 + t16 * s25 * (t306 - s13) + t392) + t185 * t226)) + t499
     & * t135 * ((t231 - t237) * s16 - t91 - t233) * t125 + t27 * (-t188
     & + t56) * t297) - t45 / ddreal(3) - ddreal(2) / ddreal(3) * t64 + 
     &t339 * t8 * t27 * t297) * struc1PM + t11 * struc7PM + ddreal(2) / 
     &ddreal(3) * t3

           ampNonresonantLightFullReC4PM = result
       end function ampNonresonantLightFullReC4PM

       function ampNonresonantLightFullReC4PP()
           implicit none
           type(dd_complex) :: ampNonresonantLightFullReC4PP

           type(dd_complex) :: result

      result = ddreal(0)

           ampNonresonantLightFullReC4PP = result
       end function ampNonresonantLightFullReC4PP

       function ampNonresonantLightFullReC7MM()
           implicit none
           type(dd_complex) :: ampNonresonantLightFullReC7MM

           type(dd_complex) :: result

      result = ddreal(0)

           ampNonresonantLightFullReC7MM = result
       end function ampNonresonantLightFullReC7MM

       function ampNonresonantLightFullReC7MP()
           implicit none
           type(dd_complex) :: ampNonresonantLightFullReC7MP

           type(dd_complex) :: result

      result = ddreal(0)

           ampNonresonantLightFullReC7MP = result
       end function ampNonresonantLightFullReC7MP

       function ampNonresonantLightFullReC7PM()
           implicit none
           type(dd_complex) :: ampNonresonantLightFullReC7PM

           type(dd_complex) :: result

      result = ddreal(0)

           ampNonresonantLightFullReC7PM = result
       end function ampNonresonantLightFullReC7PM

       function ampNonresonantLightFullReC7PP()
           implicit none
           type(dd_complex) :: ampNonresonantLightFullReC7PP

           type(dd_complex) :: result

      result = ddreal(0)

           ampNonresonantLightFullReC7PP = result
       end function ampNonresonantLightFullReC7PP

       function intHL0s25s260s56x1010D0eps0()
           implicit none
           type(dd_complex) :: intHL0s25s260s56x1010D0eps0

           type(dd_complex) :: result

      result = MBL1001(0) + (s26 + s56) * MBL1011(1)

           intHL0s25s260s56x1010D0eps0 = result
       end function intHL0s25s260s56x1010D0eps0

       function intHL0s25s260s56x1010D0eps1()
           implicit none
           type(dd_complex) :: intHL0s25s260s56x1010D0eps1

           type(dd_complex) :: result

      result = ddreal(1)

           intHL0s25s260s56x1010D0eps1 = result
       end function intHL0s25s260s56x1010D0eps1

       function intHL0s25s260s56x1011D2eps0()
           implicit none
           type(dd_complex) :: intHL0s25s260s56x1011D2eps0
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = ddreal(1) / ddreal(2)
      result = t1 * ((s25 + s26 + s56) * MBL1011(1) + MBL1001(0) + ddrea
     &l(1))

           intHL0s25s260s56x1011D2eps0 = result
       end function intHL0s25s260s56x1011D2eps0

       function intHL0s25s260s56x1011D2eps1()
           implicit none
           type(dd_complex) :: intHL0s25s260s56x1011D2eps1

           type(dd_complex) :: result

      result = ddreal(1) / ddreal(2)

           intHL0s25s260s56x1011D2eps1 = result
       end function intHL0s25s260s56x1011D2eps1

       function intHL0s25s260s56x1012D2eps0()
           implicit none
           type(dd_complex) :: intHL0s25s260s56x1012D2eps0
           type(dd_complex) ::  t1,t2

           type(dd_complex) :: result

      t1 = s26 + s56
      t2 = s25 + s26 + s56
      t1 = ddreal(1) / t1
      result = ddreal(2) * t2 * MBL1011(1) * t1 + t1 * (t2 * MBL1011(0) 
     &+ MBL1001(0))

           intHL0s25s260s56x1012D2eps0 = result
       end function intHL0s25s260s56x1012D2eps0

       function intHL0s25s260s56x1012D2eps1()
           implicit none
           type(dd_complex) :: intHL0s25s260s56x1012D2eps1
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s26 + s56
      t1 = ddreal(1) / t1
      result = t1 * ((s25 + s26 + s56) * MBL1011(1) + ddreal(1))

           intHL0s25s260s56x1012D2eps1 = result
       end function intHL0s25s260s56x1012D2eps1

       function intHL0s25s260s56x1013D4eps0()
           implicit none
           type(dd_complex) :: intHL0s25s260s56x1013D4eps0
           type(dd_complex) ::  t1,t2,t3

           type(dd_complex) :: result

      t1 = s25 + s26 + s56
      t2 = s26 + s56
      t3 = ddreal(1) / (t2)**2
      result = t3 * ((ddreal(2) * s25 + ddreal(3) * t2) * MBL1001(0) / d
     &dreal(4) + ddreal(3) / ddreal(2) * (t1)**2 * MBL1011(1)) + t1 * t3
     & * (t1 * MBL1011(0) + ddreal(1)) / ddreal(2)

           intHL0s25s260s56x1013D4eps0 = result
       end function intHL0s25s260s56x1013D4eps0

       function intHL0s25s260s56x1013D4eps1()
           implicit none
           type(dd_complex) :: intHL0s25s260s56x1013D4eps1
           type(dd_complex) ::  t1,t2,t3

           type(dd_complex) :: result

      t1 = s26 + s56
      t2 = s25 + s26 + s56
      t3 = ddreal(1) / (t1)**2
      result = t3 * (s25 / ddreal(2) + ddreal(3) / ddreal(4) * t1 + (t2)
     &**2 * MBL1011(1) / ddreal(2))

           intHL0s25s260s56x1013D4eps1 = result
       end function intHL0s25s260s56x1013D4eps1

       function intHL0s25s260s56x1020D2eps0()
           implicit none
           type(dd_complex) :: intHL0s25s260s56x1020D2eps0
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = ddreal(1) / ddreal(2)
      result = t1 * ((s26 + s56) * MBL1011(1) + MBL1001(0))

           intHL0s25s260s56x1020D2eps0 = result
       end function intHL0s25s260s56x1020D2eps0

       function intHL0s25s260s56x1020D2eps1()
           implicit none
           type(dd_complex) :: intHL0s25s260s56x1020D2eps1

           type(dd_complex) :: result

      result = ddreal(1) / ddreal(2)

           intHL0s25s260s56x1020D2eps1 = result
       end function intHL0s25s260s56x1020D2eps1

       function intHL0s25s260s56x1021D2eps0()
           implicit none
           type(dd_complex) :: intHL0s25s260s56x1021D2eps0
           type(dd_complex) ::  t1,t2

           type(dd_complex) :: result

      t1 = s26 + s56
      t2 = ddreal(1) / t1
      result = -(s25 * MBL1011(0) + (ddreal(2) * s25 + t1) * MBL1011(1) 
     &+ MBL1001(0)) * t2

           intHL0s25s260s56x1021D2eps0 = result
       end function intHL0s25s260s56x1021D2eps0

       function intHL0s25s260s56x1021D2eps1()
           implicit none
           type(dd_complex) :: intHL0s25s260s56x1021D2eps1
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s26 + s56
      t1 = ddreal(1) / t1
      result = -t1 * (s25 * MBL1011(1) + ddreal(1))

           intHL0s25s260s56x1021D2eps1 = result
       end function intHL0s25s260s56x1021D2eps1

       function intHL0s25s260s56x1022D4eps0()
           implicit none
           type(dd_complex) :: intHL0s25s260s56x1022D4eps0
           type(dd_complex) ::  t1,t2,t3

           type(dd_complex) :: result

      t1 = s26 + s56
      t2 = s25 + s26 + s56
      t3 = ddreal(1) / (t1)**2
      result = -t3 * ((MBL1001(0) + ddreal(1)) * (ddreal(2) * s25 + t1) 
     &+ t2 * (ddreal(6) * s25 + t1) * MBL1011(1)) / ddreal(2) - s25 * t2
     & * MBL1011(0) * t3

           intHL0s25s260s56x1022D4eps0 = result
       end function intHL0s25s260s56x1022D4eps0

       function intHL0s25s260s56x1022D4eps1()
           implicit none
           type(dd_complex) :: intHL0s25s260s56x1022D4eps1
           type(dd_complex) ::  t1,t2

           type(dd_complex) :: result

      t1 = s26 + s56
      t2 = ddreal(1) / (t1)**2
      result = t2 * (-s25 - t1 / ddreal(2) - s25 * (s25 + s26 + s56) * M
     &BL1011(1))

           intHL0s25s260s56x1022D4eps1 = result
       end function intHL0s25s260s56x1022D4eps1

       function intHL0s25s260s56x1031D4eps0()
           implicit none
           type(dd_complex) :: intHL0s25s260s56x1031D4eps0
           type(dd_complex) ::  t1,t2,t3

           type(dd_complex) :: result

      t1 = s26 + s56
      t2 = ddreal(2)
      t3 = ddreal(1) / (t1)**2
      result = t3 * ((s25 * t2 - t1) * MBL1001(0) - ((s26)**2 + (s56)**2
     & - t2 * ((s25 - s26) * s56 + s25 * s26) - ddreal(6) * (s25)**2) * 
     &MBL1011(1)) / ddreal(4) + s25 * t3 * (s25 * MBL1011(0) + ddreal(1)
     &) / ddreal(2)

           intHL0s25s260s56x1031D4eps0 = result
       end function intHL0s25s260s56x1031D4eps0

       function intHL0s25s260s56x1031D4eps1()
           implicit none
           type(dd_complex) :: intHL0s25s260s56x1031D4eps1
           type(dd_complex) ::  t1,t2,t3,t4

           type(dd_complex) :: result

      t1 = s26 + s56
      t2 = ddreal(2)
      t3 = ddreal(1) / ddreal(4)
      t4 = ddreal(1) / (t1)**2
      result = t4 * (t3 * (s25 * t2 - t1) + (s25)**2 * MBL1011(1) / ddre
     &al(2))

           intHL0s25s260s56x1031D4eps1 = result
       end function intHL0s25s260s56x1031D4eps1

       function intHL0s25s26s34s56x1110D2eps0()
           implicit none
           type(dd_complex) :: intHL0s25s26s34s56x1110D2eps0
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = ddreal(1) / ddreal(2)
      result = t1 * (s34 * MBL1110(1) + (s26 + s56) * MBL1011(1) + MBL10
     &01(0) + ddreal(1))

           intHL0s25s26s34s56x1110D2eps0 = result
       end function intHL0s25s26s34s56x1110D2eps0

       function intHL0s25s26s34s56x1110D2eps1()
           implicit none
           type(dd_complex) :: intHL0s25s26s34s56x1110D2eps1

           type(dd_complex) :: result

      result = ddreal(1) / ddreal(2)

           intHL0s25s26s34s56x1110D2eps1 = result
       end function intHL0s25s26s34s56x1110D2eps1

       function intHL0s25s26s34s56x1120D2eps0()
           implicit none
           type(dd_complex) :: intHL0s25s26s34s56x1120D2eps0
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s25 + s26 - s34 + s56
      t1 = ddreal(1) / t1
      result = -ddreal(2) * s34 * MBL1110(1) * t1 - (s34 * MBL1110(0) + 
     &(s26 + s56) * MBL1011(1) + MBL1001(0)) * t1

           intHL0s25s26s34s56x1120D2eps0 = result
       end function intHL0s25s26s34s56x1120D2eps0

       function intHL0s25s26s34s56x1120D2eps1()
           implicit none
           type(dd_complex) :: intHL0s25s26s34s56x1120D2eps1
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s25 + s26 - s34 + s56
      t1 = ddreal(1) / t1
      result = -t1 * (s34 * MBL1110(1) + ddreal(1))

           intHL0s25s26s34s56x1120D2eps1 = result
       end function intHL0s25s26s34s56x1120D2eps1

       function intHL0s25s26s34s56x1130D4eps0()
           implicit none
           type(dd_complex) :: intHL0s25s26s34s56x1130D4eps0
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s34 - s25 - s26 - s56
      t1 = ddreal(1) / (t1)**2
      result = -(s25 + s26 - ddreal(3) * s34 + s56) * t1 * ((s26 + s56) 
     &* MBL1011(1) + MBL1001(0)) / ddreal(4) + s34 * t1 * (s34 * MBL1110
     &(0) + ddreal(1)) / ddreal(2) + ddreal(3) / ddreal(2) * (s34)**2 * 
     &MBL1110(1) * t1

           intHL0s25s26s34s56x1130D4eps0 = result
       end function intHL0s25s26s34s56x1130D4eps0

       function intHL0s25s26s34s56x1130D4eps1()
           implicit none
           type(dd_complex) :: intHL0s25s26s34s56x1130D4eps1
           type(dd_complex) ::  t1,t2

           type(dd_complex) :: result

      t1 = s25 + s26 - s34 + s56
      t2 = -ddreal(1) / ddreal(4)
      t1 = ddreal(1) / (t1)**2
      result = t1 * (t2 * (s25 + s26 - ddreal(3) * s34 + s56) + (s34)**2
     & * MBL1110(1) / ddreal(2))

           intHL0s25s26s34s56x1130D4eps1 = result
       end function intHL0s25s26s34s56x1130D4eps1

       function intHL0s25s26s34s56x1210D2eps0()
           implicit none
           type(dd_complex) :: intHL0s25s26s34s56x1210D2eps0
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s25 + s26 - s34 + s56
      t1 = ddreal(1) / t1
      result = ((s26 + s56) * MBL1011(1) + (s25 + s26 + s56) * MBL1110(0
     &) + (s25 + s26 + s34 + s56) * MBL1110(1) + MBL1001(0)) * t1

           intHL0s25s26s34s56x1210D2eps0 = result
       end function intHL0s25s26s34s56x1210D2eps0

       function intHL0s25s26s34s56x1210D2eps1()
           implicit none
           type(dd_complex) :: intHL0s25s26s34s56x1210D2eps1
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s25 + s26 - s34 + s56
      t1 = ddreal(1) / t1
      result = t1 * ((s25 + s26 + s56) * MBL1110(1) + ddreal(1))

           intHL0s25s26s34s56x1210D2eps1 = result
       end function intHL0s25s26s34s56x1210D2eps1

       function intHL0s25s26s34s56x1220D4eps0()
           implicit none
           type(dd_complex) :: intHL0s25s26s34s56x1220D4eps0
           type(dd_complex) ::  t1,t2

           type(dd_complex) :: result

      t1 = s25 + s26 - s34 + s56
      t2 = s25 + s26 + s56
      t1 = ddreal(1) / (t1)**2
      result = -t1 * (((s26 + s56) * MBL1011(1) + MBL1001(0) + ddreal(1)
     &) * (s25 + s26 + s34 + s56) + s34 * (s34 + ddreal(5) * t2) * MBL11
     &10(1)) / ddreal(2) - s34 * t2 * MBL1110(0) * t1

           intHL0s25s26s34s56x1220D4eps0 = result
       end function intHL0s25s26s34s56x1220D4eps0

       function intHL0s25s26s34s56x1220D4eps1()
           implicit none
           type(dd_complex) :: intHL0s25s26s34s56x1220D4eps1
           type(dd_complex) ::  t1,t2

           type(dd_complex) :: result

      t1 = s25 + s26 - s34 + s56
      t2 = -ddreal(1) / ddreal(2)
      t1 = ddreal(1) / (t1)**2
      result = t1 * (t2 * (s25 + s26 + s34 + s56) - s34 * (s25 + s26 + s
     &56) * MBL1110(1))

           intHL0s25s26s34s56x1220D4eps1 = result
       end function intHL0s25s26s34s56x1220D4eps1

       function intHL0s25s26s34s56x1310D4eps0()
           implicit none
           type(dd_complex) :: intHL0s25s26s34s56x1310D4eps0
           type(dd_complex) ::  t1,t2,t3,t4

           type(dd_complex) :: result

      t1 = s25 + s26 + s56
      t2 = s25 + s26 - s34 + s56
      t3 = ddreal(3)
      t4 = ddreal(6)
      t2 = ddreal(1) / (t2)**2
      result = t2 * (((s26 + s56) * MBL1011(1) + MBL1001(0)) * (t1 * t3 
     &- s34) + (t3 * ((s25)**2 + (s26)**2 + (s56)**2) + t4 * ((s25 + s26
     &) * s56 + s25 * s26) + ddreal(4) * t1 * s34 - (s34)**2) * MBL1110(
     &1)) / ddreal(4) + t1 * t2 * (t1 * MBL1110(0) + ddreal(1)) / ddreal
     &(2)

           intHL0s25s26s34s56x1310D4eps0 = result
       end function intHL0s25s26s34s56x1310D4eps0

       function intHL0s25s26s34s56x1310D4eps1()
           implicit none
           type(dd_complex) :: intHL0s25s26s34s56x1310D4eps1
           type(dd_complex) ::  t1,t2,t3,t4

           type(dd_complex) :: result

      t1 = s25 + s26 + s56
      t2 = ddreal(3)
      t3 = s25 + s26 - s34 + s56
      t4 = ddreal(1) / ddreal(4)
      t3 = ddreal(1) / (t3)**2
      result = t3 * (t4 * (t1 * t2 - s34) + (t1)**2 * MBL1110(1) / ddrea
     &l(2))

           intHL0s25s26s34s56x1310D4eps1 = result
       end function intHL0s25s26s34s56x1310D4eps1

       function intHLs160000x0111D0eps0()
           implicit none
           type(dd_complex) :: intHLs160000x0111D0eps0

           type(dd_complex) :: result

      result = I300s16(0)

           intHLs160000x0111D0eps0 = result
       end function intHLs160000x0111D0eps0

       function intHLs160000x0111D0eps1()
           implicit none
           type(dd_complex) :: intHLs160000x0111D0eps1

           type(dd_complex) :: result

      result = I300s16(1)

           intHLs160000x0111D0eps1 = result
       end function intHLs160000x0111D0eps1

       function intHLs160000x0111D2eps0()
           implicit none
           type(dd_complex) :: intHLs160000x0111D2eps0

           type(dd_complex) :: result

      result = ddreal(3) / ddreal(2) + s16 * I300s16(1) / ddreal(2)

           intHLs160000x0111D2eps0 = result
       end function intHLs160000x0111D2eps0

       function intHLs160000x0111D2eps1()
           implicit none
           type(dd_complex) :: intHLs160000x0111D2eps1

           type(dd_complex) :: result

      result = ddreal(1) / ddreal(2)

           intHLs160000x0111D2eps1 = result
       end function intHLs160000x0111D2eps1

       function intHLs160000x0112D2eps0()
           implicit none
           type(dd_complex) :: intHLs160000x0112D2eps0

           type(dd_complex) :: result

      result = -ddreal(2) / s16 - I300s16(1)

           intHLs160000x0112D2eps0 = result
       end function intHLs160000x0112D2eps0

       function intHLs160000x0112D2eps1()
           implicit none
           type(dd_complex) :: intHLs160000x0112D2eps1

           type(dd_complex) :: result

      result = -ddreal(1) / s16

           intHLs160000x0112D2eps1 = result
       end function intHLs160000x0112D2eps1

       function intHLs160000x0113D4eps0()
           implicit none
           type(dd_complex) :: intHLs160000x0113D4eps0

           type(dd_complex) :: result

      result = -ddreal(1) / s16 / ddreal(2) - I300s16(1) / ddreal(4)

           intHLs160000x0113D4eps0 = result
       end function intHLs160000x0113D4eps0

       function intHLs160000x0113D4eps1()
           implicit none
           type(dd_complex) :: intHLs160000x0113D4eps1

           type(dd_complex) :: result

      result = -ddreal(1) / s16 / ddreal(4)

           intHLs160000x0113D4eps1 = result
       end function intHLs160000x0113D4eps1

       function intHLs160000x0121D2eps0()
           implicit none
           type(dd_complex) :: intHLs160000x0121D2eps0

           type(dd_complex) :: result

      result = ddreal(4) / s16 + I300s16(0) + ddreal(2) * I300s16(1)

           intHLs160000x0121D2eps0 = result
       end function intHLs160000x0121D2eps0

       function intHLs160000x0121D2eps1()
           implicit none
           type(dd_complex) :: intHLs160000x0121D2eps1

           type(dd_complex) :: result

      result = ddreal(2) / s16 + I300s16(1)

           intHLs160000x0121D2eps1 = result
       end function intHLs160000x0121D2eps1

       function intHLs160000x0122D4eps0()
           implicit none
           type(dd_complex) :: intHLs160000x0122D4eps0

           type(dd_complex) :: result

      result = -ddreal(3) / ddreal(2) / s16 - I300s16(1) / ddreal(2)

           intHLs160000x0122D4eps0 = result
       end function intHLs160000x0122D4eps0

       function intHLs160000x0122D4eps1()
           implicit none
           type(dd_complex) :: intHLs160000x0122D4eps1

           type(dd_complex) :: result

      result = -ddreal(1) / s16 / ddreal(2)

           intHLs160000x0122D4eps1 = result
       end function intHLs160000x0122D4eps1

       function intHLs160000x0211D2eps0()
           implicit none
           type(dd_complex) :: intHLs160000x0211D2eps0

           type(dd_complex) :: result

      result = -ddreal(2) / s16 - I300s16(1)

           intHLs160000x0211D2eps0 = result
       end function intHLs160000x0211D2eps0

       function intHLs160000x0211D2eps1()
           implicit none
           type(dd_complex) :: intHLs160000x0211D2eps1

           type(dd_complex) :: result

      result = -ddreal(1) / s16

           intHLs160000x0211D2eps1 = result
       end function intHLs160000x0211D2eps1

       function intHLs16s25s26s34s56x1111D2eps0()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1111D2eps0
           type(dd_complex) ::  t1,t2,t3,t4,t5,t6,t7,t8

           type(dd_complex) :: result

      t1 = s25 + s26 + s56
      t2 = s13 + s14 + s34
      t3 = ddreal(2)
      t4 = s16 * s25 + (-s34 + s16 + s25 + s26) * s26 + (s26 * t3 + s16 
     &+ s25 - s34 + s56) * s56
      t5 = s26 + s56
      t6 = s25 - s16
      t7 = s25 + s26 - s34 + s56
      t8 = ddreal(1) / t4
      t4 = ddreal(1) / t4
      t4 = t4 * t1
      t2 = ddreal(1) / t2
      result = t8 * (t1 * ((-t1 * MBL1111(0) + I300s16(0)) * s16 + t5 * 
     &MBL1011(0) + t7 * MBL1110(0)) + (-s25 * t6 - s26 * t6 + (s25 - s26
     &) * s34 - (s25 + s34 - s16) * s56) * MBL1101(0)) / ddreal(2) + t4 
     &* ((-t1 * MBL1111(1) + I300s16(1)) * s16 + t5 * MBL1011(1) + t7 * 
     &MBL1110(1)) + t4 * t3 * (s13 + s14 - s25 - s26 + s34 - s56) * t2

           intHLs16s25s26s34s56x1111D2eps0 = result
       end function intHLs16s25s26s34s56x1111D2eps0

       function intHLs16s25s26s34s56x1111D2eps1()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1111D2eps1
           type(dd_complex) ::  t1,t2,t3,t4

           type(dd_complex) :: result

      t1 = s25 + s26 + s56
      t2 = s13 + s14 + s34
      t3 = ddreal(2)
      t3 = s16 * s25 + (-s34 + s16 + s25 + s26) * s26 + (s26 * t3 + s16 
     &+ s25 - s34 + s56) * s56
      t4 = ddreal(1) / t3
      t3 = ddreal(1) / t3
      t2 = ddreal(1) / t2
      result = t4 * t1 * ((-t1 * MBL1111(1) + I300s16(1)) * s16 + (s26 +
     & s56) * MBL1011(1) + (s25 + s26 - s34 + s56) * MBL1110(1)) / ddrea
     &l(2) - t1 * (-s13 - s14 + s25 + s26 - s34 + s56) * t2 * t3

           intHLs16s25s26s34s56x1111D2eps1 = result
       end function intHLs16s25s26s34s56x1111D2eps1

       function intHLs16s25s26s34s56x1112D2eps0()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1112D2eps0
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s25 + s26 + s56
      t1 = ddreal(1) / t1
      result = -ddreal(1) / s16 * ((s34 * MBL1101(0) - (s34 - s25 - s26 
     &- s56) * MBL1111D2(0)) * t1 - MBL1011(0))

           intHLs16s25s26s34s56x1112D2eps0 = result
       end function intHLs16s25s26s34s56x1112D2eps0

       function intHLs16s25s26s34s56x1112D2eps1()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1112D2eps1

           type(dd_complex) :: result

      result = MBL1011(1) / s16

           intHLs16s25s26s34s56x1112D2eps1 = result
       end function intHLs16s25s26s34s56x1112D2eps1

       function intHLs16s25s26s34s56x1112D4eps0()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1112D4eps0
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t2,t3,t4
           type(dd_complex) ::  t5,t6,t7,t8,t9

           type(dd_complex) :: result

      t1 = s16 - s25
      t2 = s16 + s25
      t3 = (s25)**2
      t4 = (s16)**2
      t5 = ddreal(2) * s16
      t6 = (s34)**2
      t7 = t4 + t3 + t6
      t8 = (t1)**2
      t9 = t8 * s26
      t10 = -t5 * t3 + t7 * s25 - (s26 * t2 + ddreal(2) * t3) * s34 + (-
     &s16 * s34 - s25 * s34 - t5 * s25 + t3 + t4) * s56 + t9
      t11 = s16 * s25
      t7 = -ddreal(2) * s34 * t2 - ddreal(2) * t11 + t7
      t11 = (-s34 + t2 + s26) * s26 + (s16 + s25 + ddreal(2) * s26 - s34
     & + s56) * s56 + t11
      t12 = ddreal(3) * s25
      t4 = (t4 + t3) * s25
      t13 = s16 * t3
      t14 = t1 * s25
      t1 = t1 * s26
      t15 = s25 + s26 - s34 + s56
      t16 = ddreal(1) / t11
      t11 = ddreal(1) / t11
      t7 = ddreal(1) / t7
      t17 = ddreal(1) / ddreal(2)
      result = t17 * ((s25 + s26 + s56) * t15 * MBL1111D2(0) * t16 + ((-
     &s25 * ((s16 + ddreal(2) * s25 + s26 - s34) * s34 + (s16 - s25 + s3
     &4) * s56 + t14 + t1) * MBL1101(0) - t15 * ((s25 - s26) * s34 + (s1
     &6 - s25 - s34) * s56 + t1 + t14) * MBL1110(1)) * s34 + t10 * (-s16
     & * I300s16(1) + MBL1001(0)) + (-ddreal(2) * s25 * (t13 - t9) + t8 
     &* (s26)**2 + (t6 * (s25 + ddreal(3) * s26) + t4) * s25 + ((-(ddrea
     &l(3) * s34 + t5) * s25 + (-s34 + s16) * s16 + t3) * s56 + (t12 * s
     &34 - (ddreal(3) * s16 + ddreal(5) * s25 + ddreal(6) * s26) * s25 -
     & t5 * s26) * s34 + ddreal(2) * t4 + ddreal(2) * t9 - ddreal(4) * t
     &13) * s56 - (s25 + s26) * ((s16 + t12) * s26 + ddreal(2) * s25 * t
     &2) * s34) * MBL1011(1)) * t7 * t11) - t10 * t7 * t11

           intHLs16s25s26s34s56x1112D4eps0 = result
       end function intHLs16s25s26s34s56x1112D4eps0

       function intHLs16s25s26s34s56x1112D4eps1()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1112D4eps1

           type(dd_complex) :: result

      result = ddreal(0)

           intHLs16s25s26s34s56x1112D4eps1 = result
       end function intHLs16s25s26s34s56x1112D4eps1

       function intHLs16s25s26s34s56x1113D4eps0()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1113D4eps0
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           type(dd_complex) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           type(dd_complex) ::  t31,t32,t33,t34,t35,t36,t37,t38,t4,t5,t6,t7
           type(dd_complex) ::  t8,t9

           type(dd_complex) :: result

      t1 = s16 - s25
      t2 = ddreal(2)
      t3 = ddreal(3) * s25
      t4 = t2 * s26
      t5 = s16 + s26
      t6 = (s34)**2
      t7 = (s25)**2
      t8 = s25 * t7
      t9 = (s16)**2
      t10 = s16 * t9
      t11 = t9 + t6 + t7
      t12 = t9 + t7
      t13 = (t1)**2
      t14 = t13 * s26
      t15 = t12 * s25
      t16 = t2 * s16
      t17 = t16 * t7
      t18 = -t17 + (s25 * (s16 - t4 - t3) + (-s34 + t3 + t5) * s34) * s3
     &4 + (-t2 * s25 * (s16 + s34) + t11) * s56 + t14 + t15
      t19 = s16 + s25
      t20 = s16 * s25
      t11 = -t2 * (s34 * t19 + t20) + t11
      t21 = (-s34 + t19 + s26) * s26 + (t4 + s16 + s25 - s34 + s56) * s5
     &6 + t20
      t22 = (s26)**2
      t23 = ddreal(4) * s26
      t24 = s16 * t7
      t25 = s26 + s56
      t26 = s16 * t25
      t27 = ddreal(4) * s16
      t28 = s25 + s26
      t29 = ddreal(5) * t9
      t30 = ddreal(3) * s16
      t31 = s25 + t30
      t32 = t2 * s25
      t30 = t32 + t30
      t33 = ddreal(7) * s16
      t34 = ddreal(4) * s25
      t35 = t1 * s25
      t36 = s16 + t32
      t1 = t1 * s26
      t37 = s25 + s26 - s34 + s56
      t11 = ddreal(1) / t11
      t26 = ddreal(1) / t26
      t21 = ddreal(1) / t21
      t38 = ddreal(1) / s16
      t25 = ddreal(1) / t25
      t8 = (t2 * (t13 * t19 * t22 + t24 * t12) - ((-(-(ddreal(3) * s34 +
     & t32) * s16 + (-s34 + s25) * s25 + t9) * s56 + t2 * ((t35 + t9) * 
     &s25 - t10) + s34 * (-s34 * t30 + s16 * (ddreal(5) * s16 + ddreal(9
     &) * s26) + (ddreal(3) * s26 + t33 + t34) * s25) - ddreal(3) * t14)
     & * s56 - ddreal(4) * s25 * t10 - ddreal(4) * t14 * t19 - ((-s25 * 
     &s34 + ddreal(6) * s16 * t28 + (t23 + t3) * s25) * s34 - ddreal(3) 
     &* t22 * t31 - ddreal(3) * t8 - t4 * ((t33 + t34) * s25 + t29) - t2
     &0 * (ddreal(9) * s16 + ddreal(8) * s25)) * s34 - (t7)**2 - ddreal(
     &3) * t13 * t22 - t24 * (t32 - t33)) * s56 + t6 * ((s26 * t30 + t3 
     &* (s25 + t16)) * s26 + t17) - s25 * s26 * s34 * t6 + t14 * (s25 * 
     &(s25 + t27) + t22) - t28 * (ddreal(4) * t20 * t19 + s26 * (t31 * s
     &26 + (t27 + t3) * s25 + t29)) * s34 - ddreal(4) * t9 * t8) * MBL10
     &11(1) * t38
      t10 = t11 * t38
      t12 = ddreal(1) / ddreal(2)
      result = t12 * (t21 * (t10 * (((t1 + t35) * s25 + (s16 * t19 + s26
     & * t36 + (-s26 - t16 - t3 + s34) * s34 + ddreal(3) * t7) * s34 - (
     &(-t36 + s34) * s34 - t35) * s56) * MBL1101(0) - t37 * ((-t32 - t5 
     &+ s34) * s34 - (s16 - s25 + s34) * s56 - t1 - t35) * MBL1110(1)) *
     & s34 + t11 * ((-t2 * (s34 * (s56)**2 + t24) + ((s25 + t4) * s34 - 
     &t2 * (s25 * t19 + t22) - ddreal(3) * s26 * t19) * s34 + (t2 * (-t2
     &0 + t6) - (ddreal(3) * t19 + t23) * s34 + t7 + t9) * s56 + t14 + t
     &15) * MBL1001(0) + t8) * t25 - (t37)**2 * MBL1111D2(0) * t38 + t18
     & * I300s16(1) * t11) + (s25 + s26 + s56) * MBL1011(0) * t26) + t10
     & * t18 * t21

           intHLs16s25s26s34s56x1113D4eps0 = result
       end function intHLs16s25s26s34s56x1113D4eps0

       function intHLs16s25s26s34s56x1113D4eps1()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1113D4eps1
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s16 * (s26 + s56)
      t1 = ddreal(1) / t1
      result = t1 * ((s25 + s26 + s56) * MBL1011(1) + ddreal(1)) / ddrea
     &l(2)

           intHLs16s25s26s34s56x1113D4eps1 = result
       end function intHLs16s25s26s34s56x1113D4eps1

       function intHLs16s25s26s34s56x1114D6eps0()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1114D6eps0
           type(dd_complex) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           type(dd_complex) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           type(dd_complex) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           type(dd_complex) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           type(dd_complex) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           type(dd_complex) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           type(dd_complex) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           type(dd_complex) ::  t175,t176,t177,t178,t179,t18,t180,t181,t182,t183,t184,t185
           type(dd_complex) ::  t186,t187,t188,t189,t19,t190,t191,t192,t193,t194,t195,t196
           type(dd_complex) ::  t197,t198,t199,t2,t20,t200,t201,t202,t203,t204,t205,t206
           type(dd_complex) ::  t207,t208,t209,t21,t210,t211,t212,t213,t214,t215,t216,t217
           type(dd_complex) ::  t218,t219,t22,t220,t221,t23,t24,t25,t26,t27,t28,t29
           type(dd_complex) ::  t3,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t4
           type(dd_complex) ::  t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t5,t50
           type(dd_complex) ::  t51,t52,t53,t54,t55,t56,t57,t58,t59,t6,t60,t61
           type(dd_complex) ::  t62,t63,t64,t65,t66,t67,t68,t69,t7,t70,t71,t72
           type(dd_complex) ::  t73,t74,t75,t76,t77,t78,t79,t8,t80,t81,t82,t83
           type(dd_complex) ::  t84,t85,t86,t87,t88,t89,t9,t90,t91,t92,t93,t94
           type(dd_complex) ::  t95,t96,t97,t98,t99

           type(dd_complex) :: result

      t1 = s16 - s25
      t2 = ddreal(2) * s16
      t3 = s16 + s25
      t4 = ddreal(6) * s16
      t5 = ddreal(5) * s25
      t6 = t5 + t4
      t7 = (s16)**2
      t8 = (t7)**2
      t9 = t7 * t8
      t10 = s16 * t8
      t11 = s16 * t7
      t12 = ddreal(17) * s25
      t13 = ddreal(13) * s16
      t14 = (t12 + t13) * s25 + ddreal(2) * t7
      t15 = ddreal(19) * s25
      t16 = ddreal(5) * s16
      t17 = t16 + t15
      t18 = s25 + s26
      t19 = ddreal(15) * s25
      t20 = ddreal(22) * s16
      t21 = ddreal(19) * t7
      t22 = ddreal(30) * s25
      t23 = ddreal(32) * s16
      t24 = ddreal(11) * t11
      t25 = ddreal(31) * t7
      t26 = (ddreal(61) * s25 + t23) * s25
      t27 = ddreal(23) * s16
      t28 = ddreal(27) * s25
      t29 = t28 + t27
      t30 = s16 * s25
      t31 = (s25)**2
      t32 = (t31)**2
      t33 = (t32)**2
      t34 = s25 * t31
      t35 = t34 * t32
      t36 = t31 * t32
      t37 = s25 * t32
      t38 = t7 * t31
      t39 = ddreal(8)
      t40 = t31 + t7
      t41 = s25 * t40
      t42 = t39 * t7
      t43 = ddreal(75) * s25
      t44 = ddreal(38) * s16
      t45 = ddreal(71) * t7
      t46 = (t44 + t43) * s25 + t45
      t47 = ddreal(60) * t11 * t3 + t31 * t46
      t48 = ddreal(3)
      t49 = ddreal(200)
      t50 = t48 * s16
      t51 = ddreal(24) * t8
      t52 = (((s25 * t49 - t50) * s25 + ddreal(58) * t7) * s25 + ddreal(
     &89) * t11) * s25
      t53 = t52 + t51
      t54 = ddreal(53) * t11
      t55 = ((s16 + ddreal(162) * s25) * s25 + ddreal(24) * t7) * s25
      t56 = t55 + t54
      t57 = (ddreal(39) * s25 + t4) * s25 + t21
      t58 = (s26)**2
      t59 = (t58)**2
      t60 = s26 * t59
      t61 = s26 * t58
      t62 = t57 * s26
      t63 = t7 * t34
      t64 = ddreal(12)
      t65 = ddreal(6) * s26
      t66 = s16 * t31
      t67 = t39 * s16
      t68 = ddreal(7) * t7
      t69 = ddreal(100) * s25
      t70 = ((ddreal(102) * s16 + t69) * s25 + ddreal(113) * t7) * s25 +
     & ddreal(67) * t11
      t71 = ddreal(107) * t7
      t72 = ddreal(32) * t11
      t73 = ((ddreal(149) * s16 + ddreal(230) * s25) * s25 + t71) * s25
      t74 = t73 + t72
      t75 = ddreal(81) * s16
      t76 = ddreal(148) * s25
      t77 = (t75 + t76) * s25 + ddreal(45) * t7
      t78 = ddreal(9) * s16
      t79 = ddreal(25) * s25
      t80 = -t79 - t78
      t81 = s25 * t70
      t82 = t63 * t3
      t43 = (ddreal(83) * s16 + t43) * s25 + ddreal(49) * t7
      t83 = ddreal(32) * t7
      t84 = (ddreal(110) * s16 + ddreal(145) * s25) * s25
      t85 = t84 + t83
      t86 = ddreal(67) * s25
      t87 = ddreal(30) * s16
      t88 = t87 + t86
      t89 = s25 * t43
      t90 = t48 * s26
      t91 = t90 * t66
      t92 = t22 + t27
      t93 = ddreal(21) * s16
      t94 = ddreal(47) * s25
      t95 = t94 + t93
      t96 = s25 * t92
      t97 = t66 * t48
      t98 = t48 * s25
      t99 = ddreal(10) * s25
      t100 = ddreal(64) * s16
      t101 = ddreal(72) * t7
      t102 = ((ddreal(91) * s25 + t100) * s25 + t101) * s25
      t103 = ddreal(55) * s16
      t104 = ddreal(88) * s25
      t105 = (t103 + t104) * s25 + t25
      t106 = ddreal(4)
      t107 = ddreal(2) * s25
      t108 = t105 * t106
      t109 = (t102 + t24) * t48
      t110 = ddreal(15) * s16
      t111 = t11 * t31
      t112 = s16 * t40
      t113 = t112 * t34
      t114 = t56 * t106
      t115 = ddreal(5) * t57
      t116 = t53 * t48
      t117 = ddreal(36) * t7
      t118 = t8 * t31
      t119 = t48 * t74
      t120 = t77 * t106
      t121 = t48 * t85
      t122 = ddreal(18)
      t123 = t38 * t122
      t124 = ddreal(30) * t59
      t125 = ddreal(48) * s26
      t126 = ddreal(9) * s26
      t127 = (s34)**2
      t128 = (t127)**2
      t129 = t127 * t128
      t130 = s34 * t128
      t131 = s34 * t127
      t132 = (t1)**2
      t133 = (t132)**2
      t134 = t1 * t132
      t135 = t133 * s25
      t136 = t133 * t14
      t137 = ddreal(35) * s26
      t138 = t133 * t61
      t139 = t34 + t11
      t140 = ddreal(10) * s26 * t29
      t141 = ddreal(63) * s16
      t142 = t11 * s25
      t143 = ddreal(10) * t62
      t144 = -ddreal(10) * t80
      t145 = ddreal(60) * s26
      t146 = ddreal(72) * s26
      t147 = t122 * s26
      t148 = t8 * t34
      t149 = t133 * t17
      t150 = ddreal(5) * t10
      t151 = t106 * s26
      t152 = ddreal(6) * t129
      t153 = ddreal(31) * s16
      t154 = s16 + s26
      t155 = s16 * t154
      t156 = t32 + t8
      t157 = t30 * t40
      t158 = t132 * t18
      t12 = (-(-ddreal(28) * t157 + ddreal(6) * t128 - ((-s34 * t80 - t5
     &7) * s34 + t132 * t29) * s34 + ddreal(7) * t156 + ddreal(42) * t38
     &) * s56 + t130 * t64 - t31 * (((-ddreal(71) * s16 + t15) * s25 + d
     &dreal(94) * t7) * s25 - ddreal(46) * t11) - t150 + s25 * t8 + ((-t
     &115 * s26 + (-(ddreal(30) * t154 + t86) * s34 + (ddreal(125) * s26
     & + t75 + t76) * s25 + ddreal(45) * t155) * s34 - t54 - t55) * s34 
     &+ t132 * (s16 * (ddreal(115) * s26 + t153) + (ddreal(135) * s26 + 
     &t103 + t104) * s25)) * s34 - t137 * t133) * s56 - s25 * ((t31 * ((
     &-t103 + t12) * s25 + ddreal(52) * t7) - ddreal(23) * t8) * s25 + t
     &150) + (((((t125 + t94 + t93) * s34 - (ddreal(120) * s16 + ddreal(
     &268) * s25 + t145) * s26 - t83 - t84) * s34 + (s26 * t144 + t120) 
     &* s26 + t72 + t73) * s34 - (t143 + t114) * s26 - t51 - t52) * s34 
     &+ t132 * ((t140 + t108) * s26 + t102 + t24)) * s34 - ddreal(2) * t
     &11 * t139 - ddreal(70) * t133 * t58 - t151 * t149 - t152
      t12 = t12 * s56 - t31 * (t31 * (((t5 - t78) * s25 - t42) * s25 + d
     &dreal(22) * t11) - ddreal(13) * t10) - t48 * (t136 * s26 - t148) -
     & ddreal(6) * s25 * t9 - ddreal(6) * t149 * t58 + (((((-(t5 + t147)
     & * s34 + (ddreal(141) * s25 + t141 + t146) * s26 + t96) * s34 - ((
     &ddreal(6) * t88 + t145) * s26 + t121) * s26 - t89) * s34 + t81 + s
     &26 * (t144 * t58 + t119) + ddreal(6) * t77 * t58) * s34 - ((ddreal
     &(6) * t56 + t143) * s26 + t116) * s26 - ddreal(60) * t142 * t3 - t
     &34 * t46) * s34 + t132 * (ddreal(30) * s25 * t139 + ((ddreal(6) * 
     &t105 + t140) * s26 + t109) * s26 + t66 * (t94 + t141))) * s34 - dd
     &real(70) * t138
      t12 = t12 * s56 - t48 * (s16 * t35 + t136 * t58) - ddreal(2) * s26
     & * ((t5 + t126) * t129 + t135 * t3 * t6) + ((-t107 * t47 * s26 + (
     &t107 * t70 * s26 + (-t107 * t43 * s26 - t106 * t88 * t61 + (((t48 
     &* t95 + t125) * s26 + t107 * t92) * s26 + t97) * s34 - t110 * t34 
     &- t121 * t58 - t123 - t124) * s34 + t119 * t58 + t120 * t61 + t87 
     &* t32 - ddreal(5) * t80 * t59 + ddreal(42) * t111 + ddreal(48) * t
     &63) * s34 - t114 * t61 - t115 * t59 - t116 * t58 - t117 * t32 - dd
     &real(30) * t113 - ddreal(48) * t118) * s34 + t132 * (t110 * t32 + 
     &ddreal(27) * t111 + ddreal(30) * t63 + t107 * (t2 + t98) * ((t78 +
     & t99) * s25 + ddreal(15) * t7) * s26 + t109 * t58 + t108 * t61 + d
     &dreal(5) * t29 * t59)) * s34 - ddreal(6) * t38 * (t34 * (-s25 - s1
     &6) + t8) + ddreal(21) * t10 * t34 - t51 * t32 - t138 * (t106 * t17
     & + t137)
      t15 = s25 * t3
      t43 = t133 * t58
      t46 = t11 * t32
      t51 = t64 * t8 * t37
      t6 = -t39 * t46 * t40 - (((((s26 * (t61 * t64 + s26 * (s26 * t95 +
     & t96) + t97) - t58 * (t5 + t65) * s34) * s34 - ((s26 * t88 + t85) 
     &* s26 + t89) * t58 - ddreal(6) * t60 - ddreal(2) * t63 - t91 * t6)
     & * s34 + t39 * t82 + (((-s26 * t80 + t77) * s26 + t74) * s26 + t81
     &) * t58 + t66 * t65 * ((t5 + t67) * s25 + t68)) * s34 - t64 * t63 
     &* t40 - (s25 * t47 + ((t56 + t62) * s26 + t53) * s26) * t58 - t66 
     &* (t65 * ((ddreal(6) * t31 + t42) * s16 + ddreal(5) * t41) + t38 *
     & t39)) * s34 + t158 * (t39 * t38 * t3 + (((s26 * t29 + t25 + t26) 
     &* s26 + ((t23 + t22) * s25 + ddreal(41) * t7) * s25 + t24) * s26 +
     & t30 * ((t20 + t19) * s25 + t21)) * s26)) * s34 - t12 * s56 - t43 
     &* ((-s26 * t17 - t14) * s26 - t15 * t6) + ddreal(2) * t63 * t156 +
     & t91 * t133 * (t2 + s25) + t51 + ddreal(7) * t133 * t60
      t12 = s34 * t3
      t14 = t31 + t7 + t127 - ddreal(2) * t30 - ddreal(2) * t12
      t17 = s26 + s56
      t21 = ddreal(2) * s26
      t24 = (-s34 + t3 + s26) * s26 + (t21 + s16 + s25 - s34 + s56) * s5
     &6 + t30
      t25 = ddreal(7) * s25
      t29 = t39 * s25
      t47 = t29 + s16
      t51 = s16 * t3
      t52 = t64 * s25
      t53 = t52 + t78
      t54 = ddreal(6) * t7
      t55 = ddreal(85)
      t56 = ddreal(5) * t7
      t57 = t64 * t11
      t62 = -t31 + t7
      t70 = ddreal(72) * t34
      t72 = ddreal(13) * s16 * t62
      t73 = t72 + t70
      t74 = t31 * t48 + t7
      t75 = s25 * (((s25 * t55 - t87) * s25 + t56) * s25 + t57)
      t76 = ddreal(14) * s16
      t77 = ddreal(40) * s16
      t80 = t122 * t7
      t81 = t48 * t11
      t83 = ((t69 + t77) * s25 + t80) * s25
      t84 = ddreal(27) * s16
      t85 = ddreal(68) * s25
      t87 = t64 * t7
      t88 = (t85 + t84) * s25 + t87
      t89 = t106 * s25
      t91 = t89 + s16
      t92 = ddreal(40) * s25
      t94 = ddreal(65) * s25
      t95 = ddreal(9) * t7
      t96 = ddreal(22) * s25
      t97 = ddreal(7) * s16
      t102 = t107 + t90
      t104 = t53 / ddreal(3)
      t105 = ddreal(11) * s16
      t108 = ddreal(38) * s25
      t109 = ddreal(14) * t7
      t114 = ddreal(9) * t7 + ddreal(9) * t58
      t115 = ddreal(24) * s16
      t116 = t64 * s26
      t119 = t133 * s26
      t120 = ddreal(26) * s16
      t121 = ddreal(17) * s16
      t136 = ddreal(36) * s26
      t137 = ddreal(27) * s26
      t139 = t64 * s16
      t140 = ddreal(32) * s25
      t70 = (s25 * ((((-t153 + t29) * s25 + ddreal(44) * t7) * s25 - ddr
     &eal(26) * t11) * s25 + t106 * t8) - ddreal(6) * t130 - (((-(t140 +
     & t139 + t126) * s34 + s16 * (t139 + t126) + (t136 + t85 + t84) * s
     &25) * s34 - (ddreal(54) * t31 + t80) * s26 - t70 - t72) * s34 + t1
     &32 * (s16 * (t137 + t67) + (t108 + t136 + t121) * s25)) * s34 + t1
     &0 + t126 * t133) * s56 + ddreal(2) * s25 * t10 + ddreal(2) * t119 
     &* t47 + t129 * t48 + t34 * (((t25 - t120) * s25 + ddreal(34) * t7)
     & * s25 - ddreal(16) * t11) - (((((t96 + t116 + t78) * s34 - (ddrea
     &l(64) * s26 + t94 + t77) * s25 - t114 - t115 * s26) * s34 + (t126 
     &* t91 + ddreal(2) * t88) * s26 + t81 + t83) * s34 - (t147 * t74 + 
     &ddreal(2) * t73) * s26 - t75) * s34 + t132 * (s25 * ((t108 + t105)
     & * s25 + t109) + ((ddreal(34) * s16 + ddreal(76) * s25) * s25 + dd
     &real(16) * t7 + t126 * t104) * s26)) * s34 - t118 + ddreal(9) * t4
     &3
      t72 = t10 * t31
      t84 = (s56)**2
      t85 = s16 * t37
      t47 = t48 * ((-t106 * t157 - (t104 * t132 + t127 * (t91 - s34)) * 
     &s34 + t32 + t8 + ddreal(2) * t74 * t127 + ddreal(6) * t38) * s56 *
     & t84 + t138) + t70 * s56 + t127 * (t34 * ((t22 - t76) * s25 + t42)
     & + s26 * (s26 * t73 + t75) + ddreal(6) * t74 * t61) + t128 * (s25 
     &* ((t22 + t93) * s25 + t42) + ((t106 * (t50 + t29) + t90) * s26 + 
     &(t77 + t94) * s25 + t95) * s26) + t129 * t102 - t130 * (s25 * (t52
     & + t97) + (t96 + t78 + t65) * s26) - t131 * (t31 * ((t76 + t92) * 
     &s25 + t87) + t48 * (t61 * t91 + t142) + s26 * (s26 * t88 + t81 + t
     &83)) + t43 * t47 - t158 * (s25 * ((t50 + t52) * s25 + t54) + s26 *
     & (s26 * t53 + t39 * t51 + ddreal(26) * t31)) * s34 + t135 * (t2 + 
     &t25) * s26 + t72 + ddreal(2) * t34 * ((t34 - t11) * s25 - t8) + t8
     &5 * (t67 - t25)
      t53 = ddreal(13) * s25
      t67 = t53 + t67
      t70 = ddreal(14) * s25
      t73 = t50 + t70
      t5 = t5 + t97
      t74 = ddreal(24) * s25
      t75 = t20 + t74
      t77 = t48 * t61
      t83 = ddreal(20) * s25
      t88 = ddreal(45) * t11
      t26 = (t26 + ddreal(42) * t7) * s25 + t88
      t91 = ddreal(100) * t7
      t94 = ((ddreal(73) * s16 + ddreal(126) * s25) * s25 + t91) * s25
      t108 = t94 + t88
      t136 = ddreal(57) * s25
      t141 = ddreal(48) * s16
      t142 = ddreal(47) * t7
      t143 = (t141 + t136) * s25
      t144 = t143 + t142
      t149 = s25 * s26
      t150 = (ddreal(132) * s16 + ddreal(118) * s25) * s25 + ddreal(97) 
     &* t7
      t153 = ddreal(57) * t7
      t159 = (ddreal(130) * s16 + ddreal(116) * s25) * s25
      t160 = t159 + t153
      t161 = ddreal(53) * s25
      t162 = t141 + t161
      t163 = s25 * t150
      t164 = t106 * t31
      t165 = ddreal(52) * s16
      t166 = t165 + t136
      t167 = ddreal(52) * s25
      t168 = ddreal(36) * s16
      t169 = t168 + t167
      t27 = t31 * (t83 + t27)
      t170 = ddreal(11) * s25
      t171 = ddreal(60) * s16
      t86 = (t86 + t171) * s25
      t172 = t86 + t80
      t173 = ddreal(108) * s25
      t174 = ddreal(88) * s16
      t175 = t174 + t173
      t176 = t48 * t58
      t177 = s25 * ((ddreal(21) * s25 + t20) * s25 + t109)
      t178 = ddreal(15) * t59
      t179 = ddreal(58) * t11
      t180 = t108 * t48
      t181 = t48 * t160
      t182 = ddreal(45) * t59
      t183 = ddreal(84) * t7
      t184 = t3 * s26
      t185 = ddreal(120) * t184
      t186 = ddreal(90) * s26
      t187 = ddreal(44) * s16
      t188 = t106 * t11
      t189 = ddreal(30) * s26
      t190 = ddreal(2) * t30
      t86 = -(t122 * t128 + ddreal(5) * t156 - ddreal(20) * t157 - (((dd
     &real(45) * s26 + t161 + t141) * s34 - t142 - t143 - ddreal(60) * t
     &184) * s34 + t132 * (ddreal(15) * s26 + t28 + t20)) * s34 + ddreal
     &(30) * t38) * s56 - ((((-ddreal(53) * s16 + t70) * s25 + t101) * s
     &25 - ddreal(38) * t11) * s25 + ddreal(2) * t8) * s25 - t10 * t48 +
     & ddreal(9) * t130 + (((-(t168 + t167 + t146) * s34 + (ddreal(192) 
     &* s16 + ddreal(212) * s25 + t186) * s26 + t153 + t159) * s34 - (t1
     &06 * t144 + t185) * s26 - t88 - t94) * s34 + t132 * ((t174 + t173 
     &+ t189) * s26 + t80 + t86)) * s34 - ddreal(20) * t119
      t53 = t86 * s56 - s25 * (((((t53 - t187) * s25 + ddreal(46) * t7) 
     &* s25 - t188) * s25 - ddreal(19) * t8) * s25 + t10 * t39) - t48 * 
     &t132 * (-(ddreal(10) * t61 + s26 * (t175 * s26 / ddreal(2) + t172)
     & + t177) * s34 + t132 * t73 * s26) - ((((-t137 - t170) * s34 + (dd
     &real(156) * s26 + t136 + t165) * s25 + ddreal(108) * s26 * t154) *
     & s34 - ((ddreal(6) * t162 + t186) * s26 + t181) * s26 - t163) * s3
     &4 + s25 * (((ddreal(122) * s25 + t100) * s25 + t183) * s25 + ddrea
     &l(90) * t11) + ((ddreal(6) * t144 + t185) * s26 + t180) * s26) * t
     &127 - ddreal(30) * t43
      t53 = t53 * s56 + t106 * t37 * t62 + (((t107 * t150 * s26 + t106 *
     & t162 * t61 + (((t96 + t137) * s26 + t164) * s34 - ((t64 * t169 / 
     &ddreal(4) + t146) * s26 + t107 * t166) * s26 - t27) * s34 + t171 *
     & t34 + t181 * t58 + t182 + ddreal(40) * t32 + ddreal(52) * t38) * 
     &s34 - t106 * s26 * (s25 * t26 + t144 * t58) - t31 * (((ddreal(42) 
     &* s16 + t92) * s25 + t117) * s25 + t179) - t58 * (ddreal(60) * t3 
     &* t58 + t180)) * s34 + t132 * (t168 * t34 + t176 * t172 + t175 * t
     &61 + t177 * t65 + t178 + ddreal(20) * t32 + ddreal(32) * t38)) * s
     &34 - ddreal(20) * t138 - t43 * t48 * t73 - t66 * (((-ddreal(9) * t
     &31 + ddreal(26) * t7) * s25 - ddreal(24) * t11) * s25 + ddreal(7) 
     &* t8) - t119 * t107 * t67
      t62 = s16 * t34
      t86 = t7 * t32
      t88 = t48 * s34
      t94 = ddreal(5) * t133 * t59
      t96 = s25 + s26 + s56
      t100 = (t98 + t187) * s25 + t117
      t101 = t106 * t7
      t117 = (t107 + t105) * s25 + t101
      t136 = ddreal(9) * s25
      t137 = (t20 + t136) * s25 + t48 * t7
      t141 = t75 / ddreal(2)
      t142 = ddreal(55) * s25
      t143 = ddreal(35) * t7
      t146 = t122 * s25
      t150 = ddreal(50) * t7
      t153 = ddreal(45) * s25
      t159 = ddreal(83) * t8
      t161 = ddreal(25) * t7
      t165 = (t139 + t170) * s25 + t161
      t171 = t165 * s26
      t172 = ddreal(99) * t11
      t173 = t76 + t89
      t175 = ddreal(66) * s25
      t177 = ddreal(59) * t7
      t180 = (ddreal(140) * s16 + t175) * s25 + t177
      t181 = t98 + t97
      t185 = s25 * ((ddreal(279) * s16 + t69) * s25 + ddreal(184) * t7)
      t186 = t31 * ((ddreal(212) * s16 + t153) * s25 + ddreal(186) * t7)
      t191 = t106 * t181
      t192 = t107 + t16
      t193 = ddreal(16) * s25
      t194 = ddreal(25) * s16
      t195 = t194 + t193
      t196 = ddreal(6) * s25
      t197 = t97 + t196
      t198 = s26 * t197
      t199 = s25 * t195
      t200 = ddreal(35) * t11
      t201 = ddreal(6) * t165 * s26
      t202 = ddreal(177) * t11
      t203 = ((ddreal(355) * s16 + ddreal(144) * s25) * s25 + ddreal(336
     &) * t7) * s25
      t204 = t203 + t172
      t109 = (t29 + t121) * s25 + t109
      t205 = t173 / ddreal(2)
      t206 = t204 * t58
      t207 = s25 * (((ddreal(555) * s16 + ddreal(160) * s25) * s25 + ddr
     &eal(568) * t7) * s25 + ddreal(315) * t11)
      t208 = t180 * t106
      t209 = (t18)**2
      t45 = t132 * t209 * (t58 * (t64 * t5 * s26 + (ddreal(262) * s16 + 
     &ddreal(120) * s25) * s25 + ddreal(250) * t7) + t151 * (((ddreal(77
     &) * s16 + t136) * s25 + t45) * s25 + t200) + t30 * ((ddreal(118) *
     & s16 + t142) * s25 + t71))
      t22 = s34 * (((((((t198 * t39 + ddreal(6) * t199) * s26 + t31 * t1
     &92 * t122) * s26 + ddreal(11) * t62) * s34 - t103 * t32 - t185 * t
     &176 - ddreal(20) * t181 * t59 - t186 * t21 - t208 * t61 - ddreal(6
     &8) * t63) * s34 + t106 * s26 * (t31 * (((ddreal(199) * s16 + t22) 
     &* s25 + ddreal(238) * t7) * s25 + t202) + t206) + t207 * t176 + t1
     &24 * t109 + t64 * t205 * t60 + ddreal(110) * t85 + ddreal(162) * t
     &11 * t34 + ddreal(184) * t86) * s34 - t18 * (s25 * (t21 * ((((ddre
     &al(317) * s16 + t153) * s25 + t49 * t7) * s25 + ddreal(229) * t11)
     & * s25 + ddreal(257) * t8) + t30 * (((ddreal(144) * s16 + ddreal(1
     &10) * s25) * s25 + ddreal(118) * t7) * s25 + ddreal(188) * t11)) +
     & ((t201 + ((ddreal(433) * s16 + ddreal(294) * s25) * s25 + ddreal(
     &350) * t7) * s25 + ddreal(475) * t11) * s26 + (((ddreal(839) * s16
     & + ddreal(330) * s25) * s25 + ddreal(518) * t7) * s25 + ddreal(781
     &) * t11) * s25 + ddreal(332) * t8) * t58)) * s34 + t45)
      t45 = ddreal(20) * t7
      t49 = ddreal(24) * t9
      t71 = ddreal(100) * t8
      t103 = ddreal(101) * s16
      t136 = t189 * t5
      t210 = (((ddreal(318) * s16 + ddreal(156) * s25) * s25 + ddreal(21
     &7) * t7) * s25 + ddreal(314) * t11) * s25
      t211 = t210 + t159
      t212 = ddreal(72) * s25
      t213 = ddreal(95) * t11
      t91 = ((t212 + t103) * s25 + t91) * s25
      t214 = t91 + t213
      t215 = ddreal(10) * t214
      t216 = s25 * ((((ddreal(491) * s16 + ddreal(140) * s25) * s25 + dd
     &real(306) * t7) * s25 + ddreal(413) * t11) * s25 + ddreal(282) * t
     &8)
      t181 = ddreal(40) * t181
      t217 = ((ddreal(229) * s16 + ddreal(84) * s25) * s25 + ddreal(196)
     & * t7) * s25
      t218 = ddreal(43) * s16
      t161 = (t218 + t74) * s25 + t161
      t219 = t205 * s26
      t220 = ddreal(24) * s26
      t221 = ddreal(28) * s16
      t9 = ((t127 * t165 + t156 * t48 - t157 * t64 + t123 - ddreal(2) * 
     &s34 * (t127 * t205 + t132 * t5)) * s56 + (t31 * ((-ddreal(37) * s1
     &6 + t52) * s25 + ddreal(28) * t7) - ddreal(32) * t8) * s25 + (((-d
     &dreal(6) * (t29 + t151 + t121) * s25 - ddreal(84) * t155 + t191 * 
     &s34) * s34 + t201 + t213 + t91) * s34 - ddreal(2) * t132 * (s16 * 
     &(ddreal(42) * s26 + t194) + (t218 + t74 + t189) * s25)) * s34 + t1
     &22 * (t111 + t119) + ddreal(11) * t10) * s56 + (((((-t221 + t146) 
     &* s25 - ddreal(62) * t7) * s25 + ddreal(168) * t11) * s25 - ddreal
     &(122) * t8) * s25 + ddreal(20) * t10) * s25 + ddreal(6) * t9 + (((
     &(ddreal(140) * s16 * t18 + (t175 + t145) * s25 + t177) * s34 - ddr
     &eal(30) * s26 * (t219 + t109) - t172 - t203) * s34 + (ddreal(5) * 
     &t214 + ddreal(15) * t171) * s26 + t159 + t210) * s34 - t132 * ((dd
     &real(10) * t161 + t136) * s26 + t200 + t217)) * s34 - ddreal(2) * 
     &t197 * t130 + ddreal(45) * t43 + ddreal(5) * t119 * t141
      t9 = t9 * s56 + t35 * t64 + (((((t181 * s26 + t208) * s26 + t185) 
     &* s34 - ((ddreal(60) * t109 + ddreal(40) * t219) * s26 + t106 * t2
     &04) * s26 - t207) * s34 + ((ddreal(20) * t171 + t215) * s26 + t106
     & * t211) * s26 + t216) * s34 - t132 * (s25 * (((ddreal(311) * s16 
     &+ ddreal(64) * s25) * s25 + ddreal(312) * t7) * s25 + ddreal(129) 
     &* t11) + ((ddreal(40) * t5 * s26 + ddreal(20) * t161) * s26 + t106
     & * (t217 + t200)) * s26)) * s34 + t89 * t129 - (ddreal(2) * (t220 
     &+ t194 + t193) * s25 + ddreal(2) * t221 * s26) * t130 + s25 * t49 
     &+ t122 * s16 * t36 - ddreal(30) * t72 - ddreal(108) * t148 - ddrea
     &l(168) * t7 * t37 + ddreal(252) * t46 + ddreal(60) * t138 + ddreal
     &(10) * t43 * t141 + t119 * t39 * t137
      t9 = t9 * s56 - t48 * (t130 * (((t196 + t110) * s25 + t21 * t195) 
     &* s25 + t106 * t197 * t58 - s25 * (t151 + s25) * s34) - t33) - (((
     &ddreal(398) * s16 * t32 - (ddreal(6) * t180 * t58 + t181 * t61 + t
     &185 * t90 + ddreal(45) * t32 + ddreal(186) * t38 + ddreal(212) * t
     &62) * s34 + ddreal(60) * t109 * t61 + t124 * t205 + t207 * t90 + d
     &dreal(354) * t111 + ddreal(6) * t206 + ddreal(60) * t37 + ddreal(4
     &76) * t63) * s34 - ddreal(288) * t11 * t34 - t178 * t165 - ddreal(
     &6) * t211 * t58 - t215 * t61 - t216 * t90 - ddreal(351) * t118 - d
     &dreal(45) * t36 - ddreal(372) * t85 - ddreal(272) * t86) * s34 + t
     &158 * (s25 * (((ddreal(209) * s16 + t146) * s25 + ddreal(260) * t7
     &) * s25 + t202) + ((ddreal(10) * (ddreal(65) * s16 + ddreal(33) * 
     &s25) * s25 + ddreal(10) * t150 + t136) * s26 + ddreal(2) * (t50 + 
     &s25) * ((ddreal(87) * s25 + t103) * s25 + t143)) * s26)) * s34 + t
     &182 * t133 + ddreal(10) * t138 * t141 + t43 * t64 * t137 + t66 * (
     &((((-ddreal(122) * s16 + t140) * s25 + ddreal(108) * t7) * s25 + d
     &dreal(43) * t11) * s25 - t71) * s25 + ddreal(36) * t10) + t147 * t
     &135 * t117
      t29 = t209 * (((t171 + (ddreal(77) * t31 + ddreal(95) * t7) * s16 
     &+ ddreal(50) * t41) * s26 + (((ddreal(152) * s16 + t153) * s25 + d
     &dreal(92) * t7) * s25 + ddreal(124) * t11) * s25 + t159) * t58 + d
     &dreal(36) * t38 * t40 + t30 * (t21 * (((t168 + t142) * s25 + t143)
     & * s25 + t179) + ddreal(24) * t38)) * t127
      t15 = t132 * t18 * t209 * (s26 * (s16 * ((ddreal(46) * s16 + t142)
     & * s25 + t143) + s26 * ((t146 + t187) * s25 + t150)) + ddreal(24) 
     &* t15 * t7 + ddreal(2) * t5 * t61) * s34
      t9 = (-t130 * (t31 * (t126 * t192 + ddreal(11) * t30) + ddreal(2) 
     &* t58 * (t198 + t199)) + t62 * t133 * (t170 + t115)) * s26 + (t9 *
     & s56 + t122 * t43 * (s25 * t117 + t61) - t22 + t105 * t33 - t45 * 
     &t35 + t49 * t34 - ddreal(30) * t11 * t36 - t55 * t10 * t32 + t71 *
     & t37 + t94 * t141 + t138 * t39 * t137 + t149 * (t107 * t133 * t100
     & + t152 * (t21 + s25))) * s56 + t128 * (s26 * (((s26 * t180 + t185
     &) * s26 + t186) * s26 + t62 * (ddreal(68) * s16 + t142)) + t191 * 
     &t60 + ddreal(6) * t86) - t18 * (ddreal(24) * t82 + ((t173 * s26 + 
     &(ddreal(44) * s25 + t174) * s25 + t183) * s26 + ((ddreal(267) * s1
     &6 + t69) * s25 + ddreal(252) * t7) * s25 + t172) * t61 + t149 * (t
     &151 * (((ddreal(72) * s16 + t19) * s25 + ddreal(79) * t7) * s25 + 
     &ddreal(54) * t11) + t190 * ((ddreal(80) * s16 + t142) * s25 + ddre
     &al(69) * t7))) * t131 + t29 + t133 * t141 * t60 - t15 + s25 * (t12
     &9 * (t98 + t151) + t135 * t100) * t58
      t10 = t48 * t155
      t15 = ddreal(2) * t11
      t21 = s16 * t106
      t22 = (t21 + t196) * s25 + t7
      t29 = ((ddreal(33) * s16 + t92) * s25 + t45) * s25
      t33 = (ddreal(20) * s16 + t74) * s25 + t54
      t35 = t107 + s16
      t38 = ddreal(2) * t35
      t41 = s25 * t1
      t43 = -t34 + t11
      t45 = t134 * t58
      t46 = t31 * t43 + t45
      t49 = t134 * s26
      t55 = s25 * t43
      t69 = t62 * t1
      t71 = t66 * t1
      t35 = ((t48 * t71 + ((t127 + t22) * s34 + t41 * t104) * s34 - t55 
     &- t38 * t131) * s56 - t48 * (t55 + t49) * s25 + ((t127 * (t90 + t1
     &70 + t4) + ((ddreal(26) * s25 + t13) * s25 + t54) * s25 + t15 + t9
     &0 * t22) * s34 + t41 * (s16 * (t50 + t126) + (t21 + t70 + t116) * 
     &s25)) * s34 - ddreal(2) * t131 * (t127 + (ddreal(10) * s16 + t65 +
     & t52) * s25 + t10) + ddreal(9) * t69) * s56 - t48 * t46 * s25 - ((
     &((-t127 - ddreal(22) * s25 * t154 - (t90 + t139) * s26 - ddreal(25
     &) * t31 - t54) * s34 + (t65 * t35 + ddreal(2) * t33) * s26 + t188 
     &+ t29) * s34 - s25 * (((ddreal(35) * s25 + t97) * s25 + t87) * s25
     & + ddreal(6) * t11) - (((t167 + t120) * s25 + t87) * s25 + t188) *
     & s26 - t8 - t176 * t22) * s34 - t41 * t18 * (s16 * (t126 + t4) + (
     &t193 + t116 - s16) * s25)) * s34 + ddreal(9) * s16 * t32 * t1 - t1
     &06 * (t107 + t154) * t130 - t65 * t134 * t31
      t52 = t134 * t61
      t70 = t43 * t34
      t72 = t149 * t134
      t45 = t72 + t45 + t62 * (-s16 + s25)
      t82 = s25 + s26 - s34 + s56
      t85 = -t89 + s16
      t87 = -ddreal(2) * t85
      t90 = ddreal(42) * s25
      t91 = ddreal(21) * t184
      t16 = ((((t87 - t88) * s34 + ddreal(7) * t1 * t3) * s34 + ddreal(6
     &) * t30 * t1 - ddreal(2) * t43) * s56 + t122 * t71 + ddreal(6) * t
     &128 - ddreal(6) * t49 - ddreal(6) * t55 + ((-(t74 + t97 + t126) * 
     &s34 - s16 * (t21 + t65) + (ddreal(36) * s25 + t76 + t220) * s25) *
     & s34 + t1 * ((t74 + t13) * s25 + t56 + t91)) * s34) * s56 + t122 *
     & t69 + t48 * t128 * (t50 + t196 + t151 - s34) + ((-((t125 + t23 + 
     &t90) * s25 + t114 + t76 * s26) * s34 + ((ddreal(19) * s16 + ddreal
     &(48) * s25) * s25 + t101) * s25 + ((t212 + t221) * s25 - t42 - t65
     & * t85) * s26 + t81) * s34 + t1 * t18 * ((t16 + t28) * s25 + ddrea
     &l(10) * t7 + t91)) * s34 - ddreal(6) * t46 - t72 * t64
      t28 = ddreal(1) / s16
      t14 = ddreal(1) / (t14)**2
      t30 = (t82)**2
      t24 = ddreal(1) / (t24)**2
      t17 = ddreal(1) / (t17)**2
      t3 = (-t39 * t86 * t40 - t53 * s56 + t127 * (t64 * (t3 * t60 + t11
     &3) + ddreal(2) * t149 * (s25 * (((t93 + t83) * s25 + t80) * s25 + 
     &ddreal(29) * t11) + s26 * t26) + t61 * (s26 * t144 + t108) + t42 *
     & t32) + t128 * (t122 * t59 + s26 * ((s25 * t166 + s26 * t169) * s2
     &6 + t27) + t2 * t34) - t131 * ((ddreal(9) * t59 + t164 * ((t110 + 
     &t99) * s25 + ddreal(13) * t7)) * s26 + t39 * t51 * t34 + ((s26 * t
     &162 + t160) * s26 + t163) * t58) - s26 * ((t170 + t126) * s26 + t1
     &64) * t130 + t138 * t73 - t158 * ((t77 + t89 * (s25 * t5 + t54)) *
     & s26 + t39 * t66 * t3 + t58 * (t75 * s26 + (ddreal(43) * s25 + t44
     &) * s25 + t80)) * s34 + t135 * t67 * t58 + t119 * t31 * (t97 + t89
     &) + ddreal(2) * t62 * t156 - t88 * (-t106 * t12 + t127 * t48 - t19
     &0 + t40) * s56 * (t84)**2 + t57 * t37 + t94) * MBL1001(0) * t17
      t1 = t24 * t14 * (t47 * I300s16(1) + (t17 * ((t9 + ddreal(6) * s25
     & * (t138 * t117 + t63 * (t32 + t8)) - ddreal(24) * t11 * t37 * t40
     & + ddreal(36) * t8 * t36 + t133 * t59 * (ddreal(2) * t137 + t176))
     & * MBL1011(1) + t6) - s34 * t82 * (s34 * ((((-s34 * t102 + s25 * (
     &t97 + t99) + (t78 + t65 + t146) * s26) * s34 - s25 * ((t83 + t110)
     & * s25 + t42) - t77 - s26 * ((t74 + t97) * s26 + (t23 + t90) * s25
     & + t95)) * s34 + t18 * ((ddreal(20) * t31 + t42) * s25 + (t87 * s2
     &6 - t106 * ((-t21 - t25) * s25 + t7)) * s26 + t112 * t48)) * s34 +
     & t1 * t209 * ((t99 - s16) * s25 + ddreal(7) * t184 + t56)) + t16 *
     & s56 - ddreal(2) * t52 - ddreal(2) * t70 - ddreal(6) * s25 * t45) 
     &* MBL1110(1)) * t28 + t3)
      result = t28 * (t24 * (-t96 * t82 * t30 * MBL1111D2(0) - s34 * (-t
     &48 * t31 * t45 + (-t52 - t70) * s25 + ((((t127 * t18 + s25 * ((t19
     & + t13) * s25 + t54) + ((t170 + t4 + s26) * s26 + (t20 + t79) * s2
     &5 + t54) * s26) * s34 - s25 * (((t83 + t139) * s25 + ddreal(11) * 
     &t7) * s25 + t188) - s26 * (s26 * t33 + t188 + t29) - t38 * t61) * 
     &s34 + t18 * (s25 * (((-t2 + t19) * s25 + t68) * s25 + t188) + (t22
     & * s26 + ((t83 + t78) * s25 + t56) * s25 + t15) * s26 + t8)) * s34
     & + t41 * t209 * ((-t2 + t196 + t151) * s25 + t10)) * s34 + t35 * s
     &56 - ddreal(2) * t18 * (t2 + t98 + s26) * t130) * MBL1101(0) * t14
     &) + (t96)**2 * MBL1011(0) * t17) / ddreal(6) + t1 / ddreal(12)

           intHLs16s25s26s34s56x1114D6eps0 = result
       end function intHLs16s25s26s34s56x1114D6eps0

       function intHLs16s25s26s34s56x1114D6eps1()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1114D6eps1
           type(dd_complex) ::  t1,t2,t3

           type(dd_complex) :: result

      t1 = s26 + s56
      t2 = s25 + s26 + s56
      t3 = ddreal(1) / (t1)**2
      result = ddreal(1) / s16 * t3 * (s25 / ddreal(6) + t1 / ddreal(4) 
     &+ (t2)**2 * MBL1011(1) / ddreal(6))

           intHLs16s25s26s34s56x1114D6eps1 = result
       end function intHLs16s25s26s34s56x1114D6eps1

       function intHLs16s25s26s34s56x1121D2eps0()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1121D2eps0
           type(dd_complex) ::  t1,t2,t3

           type(dd_complex) :: result

      t1 = s25 + s26 + s56
      t2 = -s16 + s25
      t3 = ddreal(1) / s16
      t1 = ddreal(1) / t1
      result = ddreal(2) * s25 * s34 * MBL1101(0) * t3 * (t1)**2 - t1 * 
     &((s25 * MBL1011(0) + s34 * MBL1110(0) + (-s25 * t2 - s26 * t2 + (s
     &25 - s26) * s34 - (-s16 + s25 + s34) * s56) * MBL1111D2(0) * t1) *
     & t3 - I300s16(0))

           intHLs16s25s26s34s56x1121D2eps0 = result
       end function intHLs16s25s26s34s56x1121D2eps0

       function intHLs16s25s26s34s56x1121D2eps1()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1121D2eps1
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s25 + s26 + s56
      t1 = ddreal(1) / t1
      result = -t1 * ((s25 * MBL1011(1) + s34 * MBL1110(1)) / s16 - I300
     &s16(1))

           intHLs16s25s26s34s56x1121D2eps1 = result
       end function intHLs16s25s26s34s56x1121D2eps1

       function intHLs16s25s26s34s56x1121D4eps0()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1121D4eps0
           type(dd_complex) ::  t1,t2,t3,t4

           type(dd_complex) :: result

      t1 = s25 + s26 + s56
      t2 = ddreal(2)
      t2 = s16 * s25 + (-s34 + s16 + s25 + s26) * s26 + (s26 * t2 + s16 
     &+ s25 - s34 + s56) * s56
      t3 = s25 - s16
      t4 = ddreal(1) / t2
      t2 = ddreal(1) / t2
      result = t1 * t2 + t4 * (t1 * (s16 * I300s16(1) - s34 * MBL1110(1)
     & - t1 * MBL1011(1) - MBL1001(0)) + (-s25 * t3 - s26 * t3 + (s25 - 
     &s26) * s34 - (s25 - s16 + s34) * s56) * MBL1111D2(0) - s25 * s34 *
     & MBL1101(0)) / ddreal(2)

           intHLs16s25s26s34s56x1121D4eps0 = result
       end function intHLs16s25s26s34s56x1121D4eps0

       function intHLs16s25s26s34s56x1121D4eps1()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1121D4eps1

           type(dd_complex) :: result

      result = ddreal(0)

           intHLs16s25s26s34s56x1121D4eps1 = result
       end function intHLs16s25s26s34s56x1121D4eps1

       function intHLs16s25s26s34s56x1122D4eps0()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1122D4eps0
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t2,t3,t4,t5,t6,t7
           type(dd_complex) ::  t8,t9

           type(dd_complex) :: result

      t1 = s25 + s26 - s34 + s56
      t2 = s16 + s25
      t3 = ddreal(2)
      t4 = s16 * s25 + (-s34 + t2 + s26) * s26 + (s26 * t3 + s16 + s25 -
     & s34 + s56) * s56
      t5 = s25 + s26 + s56
      t6 = s26 + s56
      t7 = s16 * t6
      t8 = ddreal(3) * s16
      t9 = (s56)**2
      t10 = (s26)**2
      t11 = -s16 + s25
      t12 = t11 * s26
      t7 = ddreal(1) / t7
      t4 = ddreal(1) / t4
      t6 = ddreal(1) / t6
      t13 = ddreal(1) / s16
      t14 = ddreal(1) / t5
      result = -t3 * t1 * t13 * t4 - t4 * (-s34 * t13 * (s25 * MBL1101(0
     &) * t14 + MBL1110(1)) + I300s16(1)) * t1 - t4 * (t5 * MBL1001(0) +
     & (t3 * s16 * (s25)**2 + (t3 * t2 * s26 + s25 * (s25 - s34) + t8 * 
     &s25) * s56 + t2 * t9 + t2 * t10 - s25 * (s34 - s25 - t8) * s26) * 
     &MBL1011(1) * t13) * t6 - s25 * MBL1011(0) * t7 + (t3 * ((-s34 * (s
     &25 + s26) + t12) * s25 + ((-s16 - s34 + s25) * s25 + t12) * s56) +
     & (s25 * t11 + (s34)**2) * s25 + t11 * (t10 + t9)) * MBL1111D2(0) *
     & t13 * t14 * t4

           intHLs16s25s26s34s56x1122D4eps0 = result
       end function intHLs16s25s26s34s56x1122D4eps0

       function intHLs16s25s26s34s56x1122D4eps1()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1122D4eps1
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s16 * (s26 + s56)
      t1 = ddreal(1) / t1
      result = -t1 * (s25 * MBL1011(1) + ddreal(1))

           intHLs16s25s26s34s56x1122D4eps1 = result
       end function intHLs16s25s26s34s56x1122D4eps1

       function intHLs16s25s26s34s56x1123D6eps0()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1123D6eps0
           type(dd_complex) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           type(dd_complex) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           type(dd_complex) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           type(dd_complex) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           type(dd_complex) ::  t142,t143,t144,t145,t146,t147,t148,t15,t16,t17,t18,t19
           type(dd_complex) ::  t2,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3
           type(dd_complex) ::  t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40
           type(dd_complex) ::  t41,t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51
           type(dd_complex) ::  t52,t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62
           type(dd_complex) ::  t63,t64,t65,t66,t67,t68,t69,t7,t70,t71,t72,t73
           type(dd_complex) ::  t74,t75,t76,t77,t78,t79,t8,t80,t81,t82,t83,t84
           type(dd_complex) ::  t85,t86,t87,t88,t89,t9,t90,t91,t92,t93,t94,t95
           type(dd_complex) ::  t96,t97,t98,t99

           type(dd_complex) :: result

      t1 = ddreal(2)
      t2 = t1 * s25
      t3 = s16 - s25
      t4 = ddreal(3)
      t5 = t4 * s25
      t6 = s16 - t5
      t7 = s25 + s26
      t8 = (s16)**2
      t9 = (t8)**2
      t10 = s16 * t8
      t11 = t4 * s16
      t12 = ddreal(18) * s25
      t13 = s16 + s25
      t14 = ddreal(19)
      t15 = ddreal(4)
      t16 = t14 * s16
      t17 = ddreal(20) * s25
      t18 = ddreal(7) * t8
      t19 = t13 * s26
      t20 = ddreal(21)
      t21 = ddreal(6)
      t22 = s16 * t21
      t23 = (s25 * t20 + t22) * s25
      t24 = ddreal(8)
      t25 = ddreal(17) * s25
      t26 = t24 * s16
      t27 = (s26)**2
      t28 = (t27)**2
      t29 = s26 * t27
      t30 = t8 * s25
      t31 = (s25)**2
      t32 = (t31)**2
      t33 = s25 * t32
      t34 = s25 * t31
      t35 = s16 * t31
      t36 = ddreal(30)
      t37 = s16 + s26
      t38 = s26 * t37
      t39 = ddreal(7) * s16
      t40 = t1 * s26
      t41 = ddreal(5) * s25
      t42 = s16 * t13
      t43 = ddreal(22) * t31
      t44 = ddreal(11)
      t45 = t44 * s16
      t46 = t15 * t8
      t47 = t4 * s26
      t48 = t15 * s16
      t49 = t24 * s26
      t50 = (t3)**2
      t51 = t10 * s25
      t52 = t50 * s26
      t53 = (s34)**2
      t54 = (t53)**2
      t55 = s34 * t53
      t56 = ddreal(12) * s25
      t57 = t21 * s26
      t58 = ddreal(9) * s25
      t59 = t52 * t21
      t60 = s34 * t13
      t61 = s16 * s25
      t62 = -t1 * (t60 + t61) + t8 + t31 + t53
      t63 = t10 * t31
      t64 = t50 * t27
      t65 = (s56)**2
      t66 = t50 * t29
      t67 = ddreal(13) * s16
      t68 = t31 * ((t67 - t41) * s25 - t44 * t8) + t59 * (s16 - t2)
      t69 = (-s34 + t13 + s26) * s26 + (s16 + s25 - s34 + t40 + s56) * s
     &56 + t61
      t70 = ddreal(5) * s16
      t71 = s16 + t2
      t72 = ddreal(7) * s25
      t73 = s16 + t72
      t74 = t8 + t31
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
      t86 = ddreal(36) * s16
      t87 = ddreal(9) * t8
      t88 = (ddreal(51) * s25 + t86) * s25 + t87
      t89 = s16 * t34
      t90 = ddreal(10) * t89
      t91 = t21 * s25
      t92 = t8 * t31
      t93 = t15 * s26
      t94 = s25 + t47
      t95 = s16 * t32
      t96 = ddreal(33) * s25
      t97 = ddreal(9) * s16
      t98 = -t64 + t89
      t99 = t15 * t10
      t100 = t4 * t8
      t101 = ddreal(5) * t8
      t102 = t1 * t74
      t103 = t48 * s25
      t104 = t89 * t74
      t105 = t50 * t28
      t106 = t8 * t34
      t78 = ((((s34 - t11 - t41) * s34 + t102 - t103) * s56 - s25 * ((t6
     &7 - t72) * s25 - t101) - ((-t11 - t93 - t81 + s34) * s34 + ddreal(
     &12) * s16 * t7 + (ddreal(20) * s26 + t25) * s25 + t100) * s34 + t1
     &0 + t52 * t24) * s56 - (-t24 * t34 - t99) * s25 + t4 * (-(s25 * t8
     &5 + s26 * (t1 * t78 * s26 + t88 / ddreal(3))) * s34 + t52 * t73) -
     & ddreal(12) * t98 + t53 * (-(t47 + t41) * s34 + (t96 + t57 + t97) 
     &* s26 + t83)) * s56 + t24 * t52 * (s25 * t71 + t27) + t4 * (t64 * 
     &t73 + t33) - (t91 * t85 * s26 + (-((t4 * t82 + t93) * s26 + t83 * 
     &t1) * s26 - t31 * (t58 + t45) + (s26 + t5) * t94 * s34) * s34 + t8
     &0 * t15 + t27 * t88 + ddreal(9) * t32 + t90 + ddreal(13) * t92) * 
     &s34 - t95 - t92 * (-t70 + t72)
      t85 = s34 * t7
      t88 = s26 + s56
      t93 = s16 - t77
      t107 = -s16 * t3
      t108 = t24 * s25
      t109 = (t22 + t12) * s25
      t110 = t71 * t27
      t111 = s25 * t7
      t112 = t61 * t4
      t113 = s26 + t2
      t114 = ddreal(10) * t31
      t115 = ddreal(16) * s25
      t116 = ddreal(12) * s16
      t117 = t19 * t21
      t118 = t21 * t8
      t119 = t113 * s34
      t120 = ddreal(10) * s25
      t121 = s16 + t120
      t122 = t84 * s26
      t123 = ddreal(25) * s25
      t124 = ddreal(17) * s16
      t125 = t123 + t124
      t126 = t115 + t70
      t67 = t31 * (t67 + t56)
      t127 = t15 * t31
      t128 = ddreal(29) * s25
      t129 = (t128 + t124) * s25 + ddreal(14) * t8
      t130 = (t48 + t91) * s25 + t8
      t131 = s25 * (t48 + t81)
      t132 = ddreal(48)
      t133 = ddreal(5) * s26
      t134 = t8 * t24
      t135 = t1 * t53
      t84 = (s56 * (-s56 * (s34 * (-s34 * t1 + t70 + t72) - t75 + t91 * 
     &s16) - s25 * (s25 * (-t120 + t16) - t134) - s34 * (t15 * (s16 * (s
     &16 + t133) + s25 * (ddreal(7) * s26 + t48 + t91)) + t135) + t53 * 
     &(t115 + t70 + t49) + t10 + ddreal(12) * t52) - (-t34 * t44 - t99) 
     &* s25 + t4 * (t121 * t52 + t92) - ddreal(18) * t98 - (((t57 + t72)
     & * s34 - s25 * (t123 + t124) - (s25 * t132 + ddreal(15) * s16 + dd
     &real(12) * s26) * s26) * s34 + s25 * t129 + (t122 * t21 + ddreal(1
     &2) * t130) * s26) * s34) * s56 - t1 * (t113 * t55 * t94 - t131 * t
     &52) - t31 * ((-t127 + t118) * s25 - ddreal(5) * t10) - t4 * (-t121
     & * t64 + t95) - s34 * (t2 * t129 * s26 + t84 * t15 * t29 - (((t126
     & * t4 + t49) * s26 + t2 * t125) * s26 + t67) * s34 + ddreal(12) * 
     &t130 * t27 + ddreal(12) * t32 + t90 + ddreal(14) * t92) + ddreal(1
     &2) * t66
      t90 = s25 * (s25 + t116) + t18
      t94 = ddreal(14) * t61
      t98 = t75 + t94
      t113 = ddreal(20) * s16
      t120 = s25 * (t113 + t58) + t8
      t94 = ddreal(5) * t74 + t94
      t123 = ddreal(15) * s25
      t124 = t13 * t1
      t129 = t2 + t70
      t130 = t5 + t48
      t136 = t71 * s26
      t137 = t21 * t31
      t44 = t89 * t44
      t138 = ddreal(5) * t94
      t139 = ddreal(42) * t8
      t140 = s25 * (ddreal(91) * s16 + t96) + t139
      t86 = (t86 + t25) * s25
      t141 = t101 + t86
      t142 = t5 * t140
      t143 = t141 * t15
      t144 = t15 * t120
      t128 = s25 * (s25 * (ddreal(95) * s16 + t128) + ddreal(83) * t8) +
     & ddreal(33) * t10
      t145 = t8 * t20
      t146 = s25 * (s25 * (ddreal(26) * s16 + t81) + t145) + t1 * t10
      t147 = ddreal(10) * t94
      t148 = t3 * s25
      t9 = s56 * (-s56 * (s34 * (-t124 * s34 + t94) + t4 * (s25 * (t8 + 
     &t148) - t10)) - s25 * (s25 * (t36 * t8 - ddreal(9) * t31) - ddreal
     &(18) * t10) + t1 * (-t55 * t71 + t89) + s34 * (s34 * (ddreal(10) *
     & t19 + t101 + t86) - s25 * (s25 * (ddreal(52) * s16 + ddreal(22) *
     & s25) + t139) - t99 - t138 * s26) + t9 + ddreal(15) * t19 * t50) -
     & s25 * (-t1 * t54 + ddreal(66) * t92) + ddreal(9) * s25 * (t9 + t3
     &2) + ddreal(24) * t35 * t74 - (s34 * (s34 * (s25 * (ddreal(16) * s
     &26 + t113 + t123) + t26 * s26) - s25 * t140 - s26 * (t143 + ddreal
     &(20) * t19)) + s25 * t128 + s26 * (t147 * s26 + t146 * t24)) * s34
     & + t52 * (t19 * t36 + t144)
      t9 = s56 * t9 - s25 * (-ddreal(9) * t52 * t98 + t30 * (t132 * t31 
     &- t145)) - t21 * (t10 * t34 - t64 * t120) + t36 * (s16 * t33 + t66
     & * t13) + t4 * (t31 * t32 + t55 * (-s25 * (t133 * t130 + ddreal(10
     &) * t61) - t15 * (t110 + t34) + s25 * (s25 + t40) * s34)) - s34 * 
     &(t5 * t128 * s26 - s34 * (t141 * t21 * t27 + t142 * s26 + ddreal(2
     &0) * t13 * t29 + ddreal(18) * t32 + ddreal(90) * t89 + ddreal(72) 
     &* t92) + ddreal(12) * t146 * t27 + t147 * t29 + ddreal(72) * t106 
     &+ ddreal(12) * t33 + ddreal(66) * t63 + ddreal(90) * t95)
      t22 = (t7)**2 * (ddreal(12) * t92 * t13 + s26 * (t27 * t94 + t61 *
     & (s25 * (ddreal(14) * s16 + t96) + ddreal(25) * t8)) + t15 * (s25 
     &* (s25 * (t22 + t5) + t134) + t10) * t27) * s34
      t86 = t7 * (s26 * (s26 * (s25 * (s25 * (ddreal(57) * s16 + t12) + 
     &ddreal(37) * t8) + s26 * (s25 * (ddreal(34) * s16 + t123) + t101))
     & + t35 * (ddreal(35) * s16 + t96)) + t118 * t34 + t124 * t28) * t5
     &3
      t9 = s26 * (-t55 * (s26 * (s26 * (t136 * t1 + t41 * t130) + t137 *
     & t129) + t44) + t89 * t50 * (t81 + t16)) + s56 * (s56 * t9 + t21 *
     & s25 * s26 * (s25 * t50 * t90 + t54 * t7) - (s34 * (s34 * (s26 * (
     &s26 * (t123 * t130 + t136 * t24) + ddreal(12) * t31 * t129) + t44)
     & - ddreal(36) * t19 * t31 * (s25 + t48) - t142 * t27 - ddreal(33) 
     &* t95 - ddreal(41) * t106 - ddreal(10) * t13 * t28 - t143 * t29) +
     & t7 * (s26 * (s25 * (s25 * (s25 * (ddreal(147) * s16 + ddreal(24) 
     &* s25) + ddreal(106) * t8) + ddreal(83) * t10) + s26 * (s25 * (s25
     & * (ddreal(138) * s16 + ddreal(63) * s25) + ddreal(143) * t8) + dd
     &real(16) * t10)) + t138 * t29 + t35 * (s25 * (ddreal(38) * s16 + t
     &96) + ddreal(49) * t8))) * s34 + t58 * t64 * t98 + t144 * t66 + dd
     &real(15) * t105 * t13 - t89 * (s25 * (s25 * (t11 - t81) + ddreal(2
     &7) * t8) - t10 * t14)) + s25 * t27 * (t40 + t5) * t54 - t22 + t86 
     &+ t105 * t120
      t22 = t21 * t111
      t29 = t74 * s25
      t44 = t74 * t31
      t50 = s25 * (-t52 + t35)
      t54 = s25 + s26 - s34 + s56
      t81 = t3 * s26
      t86 = ddreal(1) / s16
      t92 = ddreal(1) / t62
      t69 = ddreal(1) / (t69)**2
      t88 = ddreal(1) / (t88)**2
      t9 = (-t21 * t8 * t32 * t74 - t4 * t64 * (s26 * (s25 * t98 + t13 *
     & t27) + t31 * t90) - t9 + ddreal(12) * t10 * t33) * MBL1011(1) * t
     &86
      t9 = (-t1 * t104 - t84 * s56 - t53 * (t1 * (t89 + t28) + s26 * ((s
     &25 * t125 + s26 * t126) * s26 + t67)) + t85 * (t27 * ((t25 + t45) 
     &* s25 + t46 + t122) + t79 + t2 * ((t11 + t91) * s25 + t101) * s26)
     & + s26 * ((t40 + t72) * s26 + t127) * t55 - t52 * ((s26 * t121 + t
     &131) * s26 + t31 * (t70 + t77)) + t46 * t32 - t105 * t4) * MBL1001
     &(0) - t1 * (t104 + t105) + t15 * s25 * (-t64 * t71 + t106) - t78 *
     & s56 - t53 * (s26 * (((s26 + t82) * s26 + t83) * s26 + t31 * (t58 
     &+ t45)) + t76 * t34) + t85 * (t4 * s26 * (s25 * (s25 * t76 + t75) 
     &+ s26 * ((t11 + t77) * s25 + t8)) + t79 + t80) + s26 * ((s26 + t41
     &) * s26 + t31 * t4) * t55 - t52 * (t27 * t73 + t31 * (t5 + t70)) +
     & t9
      t5 = s34 * t54 * (t1 * (t44 + t64) - t15 * t50 - s34 * (s34 * (-s2
     &6 * (s26 + t76) - t112 - t22 + t119) + t7 * (t19 * t4 - t107 + t13
     &7)) - s56 * (-s56 * (-t4 * t60 + t102 - t103 + t53) - t15 * (t52 +
     & t29) + s34 * (s25 * (t76 + t58) + t117 + t8 + t53) + t35 * t24 - 
     &t135 * (t37 + t5))) * MBL1110(1)
      t5 = t5 + s25 * t68 - t1 * (s56 * t62 * t65 + t66) + t4 * (t6 * t6
     &4 + t63) + ((((-t40 - t41) * s34 + (ddreal(18) * s26 + t39 + t17) 
     &* s25 + t15 * t38) * s34 - t1 * ((t8 + t23 + t27) * s26 + t30) - t
     &36 * t34 - (t25 + t26) * t27 - t35) * s34 + t7 * (s25 * ((-t16 + t
     &17) * s25 + t18) + ((-t11 + t12) * s25 + t8 + t19 * t15) * s26)) *
     & s34 + ((((s16 * t20 - t58) * s25 - ddreal(15) * t8) * s25 + t10 *
     & t4 + t15 * t55 + s34 * (-(t25 + t57 + t26) * s34 + ddreal(12) * t
     &19 + t42 + t43) - t59) * s56 - t1 * (((t25 + t47 + t26) * s26 + t2
     &3 + t8 + t53) * s34 - s25 * ((s25 * t14 - t45) * s25 + t46) - (t19
     & * t21 + t42 + t43) * s26) * s34 + t21 * (t51 + t52 * (-s26 + t6))
     & + t31 * ((s16 * t36 - t56) * s25 - ddreal(24) * t8) + (t48 + t49 
     &+ t12) * t55) * s56
      t4 = t92 * t69 * ((t34 * ((-t2 + t70) * s25 - t46) + (((t1 * t38 +
     & t111 * t24 + t112 - t119) * s34 - ddreal(12) * t34 - (t8 + t109 +
     & t27) * s26 - t61 * t13 - t110 * t15) * s34 + t7 * (s25 * ((t108 -
     & t39) * s25 + t100) + (t1 * t19 + t24 * t31 - t107) * s26)) * s34 
     &+ ((-s56 * t62 + s25 * ((-t77 + t97) * s25 - t118) + t1 * t55 + (-
     &(t47 + t48 + t108) * s34 + t114 + t117 + t42) * s34 + t10 - t52 * 
     &t4) * s56 + t1 * ((t53 * (s16 + t40 + t77) + s25 * ((-t48 + t108) 
     &* s25 + t1 * t8) + (t19 * t4 + t114 + t42) * s26) * s34 + t51 + t5
     &2 * t93) + t31 * ((t116 - t41) * s25 - t87) - t53 * ((t115 + t47 +
     & t26) * s26 + t109 + t8 + t53) - t64 * t4) * s56 + t63 - t66 + t52
     & * (t76 - t41) * s25 + t64 * t93) * I300s16(1) + t86 * t5 + t9 * t
     &88)
      result = t4 / ddreal(4) - t86 * (t69 * (-s25 * s34 * (-t1 * (t55 *
     & (t37 + t2) + t50) + s34 * (s34 * (s16 * (t2 + t47) + t27 + t8 + t
     &22 + t53) - t7 * (s25 * (-t76 + t77) + t136 + t8)) + s56 * (s56 * 
     &(-t1 * (s16 + s34) * s25 - s34 * (s16 - s34) + t31 + t8) + t1 * (t
     &52 + t29 - t55) - s34 * (-s34 * (t40 + t11 + t91) + s26 * (t76 + t
     &77) - t107 + t137) - t35 * t15) + t64 + t44) * MBL1101(0) * t92 + 
     &t54 * (t1 * (s25 * (t81 + t85) + s56 * ((s16 + s34 - s25) * s25 + 
     &t81)) - s25 * (-t148 + t53) - t3 * (-t27 - t65)) * MBL1111D2(0)) +
     & s25 * (s25 + s26 + s56) * MBL1011(0) * t88) / ddreal(2)

           intHLs16s25s26s34s56x1123D6eps0 = result
       end function intHLs16s25s26s34s56x1123D6eps0

       function intHLs16s25s26s34s56x1123D6eps1()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1123D6eps1
           type(dd_complex) ::  t1,t2

           type(dd_complex) :: result

      t1 = s26 + s56
      t2 = ddreal(1) / (t1)**2
      result = ddreal(1) / s16 * t2 * (-s25 / ddreal(2) - t1 / ddreal(4)
     & - s25 * (s25 + s26 + s56) * MBL1011(1) / ddreal(2))

           intHLs16s25s26s34s56x1123D6eps1 = result
       end function intHLs16s25s26s34s56x1123D6eps1

       function intHLs16s25s26s34s56x1131D4eps0()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1131D4eps0
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           type(dd_complex) ::  t20,t21,t22,t23,t24,t25,t3,t4,t5,t6,t7,t8
           type(dd_complex) ::  t9

           type(dd_complex) :: result

      t1 = s16 + s25 - s34
      t2 = s26 + s56
      t3 = ddreal(2)
      t4 = t2 * t3 + t1
      t5 = s16 + s25
      t6 = s16 * s25
      t1 = (-s34 + t5 + s26) * s26 + (s26 * t3 + s56 + t1) * s56 + t6
      t7 = s25 + s26 + s56
      t8 = (s26)**2
      t9 = (s56)**2
      t10 = (s25 - s34) * s25
      t11 = s25 + s26
      t12 = s25 * s26
      t13 = s25 + s26 - s34 + s56
      t14 = (s34)**2
      t15 = (s25)**2
      t16 = s26 * s34
      t17 = t3 * s16
      t18 = s16 - s25
      t19 = t18 * s25
      t20 = t18 * s26
      t21 = (s16)**2
      t22 = t21 + t15
      t18 = (t18)**2
      t23 = t18 * s26
      t24 = ddreal(1) / t1
      t1 = ddreal(1) / t1
      t25 = ddreal(1) / s16
      t7 = ddreal(1) / t7
      t13 = ddreal(1) / t13
      t2 = ddreal(1) / t2
      t18 = (-t3 * ((-s16 * t15 + t23) * s25 + (-t17 * t15 + s25 * t22 +
     & (-s25 * t5 - t17 * s26 + t16) * s34 + t23) * s56 - t11 * (s16 * s
     &26 + t15) * s34) - t14 * (t8 + t15) - t15 * t22 - t8 * t18 - t9 * 
     &(-t3 * s16 * (s25 + s34) + t21 + t15 + t14)) * MBL1111D2(0) * t25 
     &* t7 * t1 + I300s16(0)
      t21 = ddreal(1) / ddreal(2)
      result = t21 * (t13 * (t14 * MBL1110(0) * t25 * t7 + s34 * ((s25 *
     & t3 + s16 + s26 - s34) * s34 + (s16 - s25 + s34) * s56 + t19 + t20
     &) * MBL1110(1) * t1 * t25 + t1 * ((t3 * (s56 * t11 + t12) + t10 + 
     &t8 + t9) * MBL1001(0) + (t17 * s25 * t15 + s16 * s26 * t8 + s16 * 
     &s56 * t9 + ((t16 - t3 * (s26 * t5 + t6) - t8) * s34 + t8 * (ddreal
     &(3) * s16 + s25)) * s25 + (s16 * ((ddreal(6) * s25 + ddreal(3) * s
     &26) * s26 + ddreal(4) * t15) - t3 * s25 * ((s16 + s25 + s26) * s34
     & - t12) + s25 * (t15 + t14)) * s56 + t9 * (ddreal(3) * s16 * t11 +
     & t10) + t15 * (ddreal(4) * s16 + s25) * s26) * MBL1011(1) * t25) *
     & t2) + t7 * t18 + s25 * s34 * ((s25 - s26) * s34 + (s16 - s25 - s3
     &4) * s56 + t19 + t20) * MBL1101(0) * t25 * (t7)**2 * t1 + t4 * I30
     &0s16(1) * t24 + t15 * MBL1011(0) * t25 * t2 * t7) + t4 * t25 * t1

           intHLs16s25s26s34s56x1131D4eps0 = result
       end function intHLs16s25s26s34s56x1131D4eps0

       function intHLs16s25s26s34s56x1131D4eps1()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1131D4eps1
           type(dd_complex) ::  t1,t2,t3,t4,t5,t6

           type(dd_complex) :: result

      t1 = s26 + s56
      t2 = ddreal(2)
      t3 = s25 + s26 - s34 + s56
      t4 = s25 + s26 + s56
      t4 = ddreal(1) / t4
      t3 = ddreal(1) / t3
      t5 = ddreal(1) / t1
      t6 = ddreal(1) / ddreal(2)
      result = t6 * (I300s16(1) * t4 + (t5 * ((t1 * t2 + s25 - s34) * t3
     & + (s25)**2 * MBL1011(1) * t4) + (s34)**2 * MBL1110(1) * t4 * t3) 
     &/ s16)

           intHLs16s25s26s34s56x1131D4eps1 = result
       end function intHLs16s25s26s34s56x1131D4eps1

       function intHLs16s25s26s34s56x1132D6eps0()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1132D6eps0
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           type(dd_complex) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           type(dd_complex) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           type(dd_complex) ::  t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52
           type(dd_complex) ::  t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63
           type(dd_complex) ::  t64,t65,t7,t8,t9

           type(dd_complex) :: result

      t1 = s25 + s26 + s56
      t2 = ddreal(3)
      t3 = ddreal(2)
      t4 = t3 * s16
      t5 = t2 * s25
      t6 = (s26)**2
      t7 = (t6)**2
      t8 = s26 * t6
      t9 = (-s25 + s34) * s25
      t10 = ddreal(4)
      t11 = t2 * s26
      t12 = t10 * s25
      t13 = t10 * s26
      t14 = t4 * s25
      t15 = s26 + s56
      t16 = t3 * s26
      t17 = s16 * s25
      t18 = (-s34 + s16 + s25 + s26) * s26 + (t16 + s25 + s16 - s34 + s5
     &6) * s56 + t17
      t19 = ddreal(8) * s16
      t20 = ddreal(5)
      t21 = ddreal(10) * s25
      t22 = s25 + s26
      t23 = t20 * s25
      t24 = s16 * t22
      t25 = -s16 + s34
      t26 = (s25)**2
      t27 = (t26)**2
      t28 = s25 * t26
      t29 = ddreal(10) * s26
      t30 = (s56)**2
      t31 = (t30)**2
      t32 = s56 * t30
      t33 = t3 * s25
      t34 = s26 + s16
      t35 = (-s26 + t33) * s34
      t36 = t4 * t26
      t37 = ddreal(9) * s16
      t38 = s25 + t37
      t39 = s16 + t5
      t40 = t10 * s16
      t41 = s25 + t40
      t42 = (s16)**2
      t43 = (t5 + t19) * s25
      t44 = t43 - t42
      t45 = ddreal(6) * s25
      t46 = s16 * s26
      t47 = s16 * t26
      t48 = t20 * s16
      t49 = t40 * s26
      t50 = (s34)**2
      t51 = ddreal(18) * t26 * t41
      t52 = -t13 * t44
      t53 = t3 * t26
      t54 = s16 * t6
      t55 = s16 * t2 + t45
      t56 = s16 * t25
      t57 = s34 * t22
      t58 = s25 - s16
      t59 = t58 * s26
      t60 = t58 * s25
      t61 = (t58)**2
      t62 = t61 * s26
      t63 = t61 * t6
      t64 = (t42 + t26) * s25
      t65 = ddreal(1) / s16
      t15 = ddreal(1) / (t15)**2
      t18 = ddreal(1) / (t18)**2
      t7 = (s16 * s26 * t7 + s16 * s56 * t31 + (s16 * t28 * (ddreal(14) 
     &* s16 + ddreal(11) * s25) + t57 * ((-s25 * (t40 + t45) + t46) * s2
     &6 - ddreal(11) * t47)) * s26 + ((t20 * t7 + ddreal(11) * t27) * s1
     &6 + s26 * ((-t52 + t51) * s26 + t53 * t38 * t39) - (((ddreal(9) * 
     &s25 * (s16 + t33) - t49) * s26 + ddreal(6) * t26 * (t48 + t33)) * 
     &s26 + ddreal(11) * s16 * t28) * s34 + ddreal(6) * t26 * s26 * t50 
     &+ ddreal(14) * t42 * t28) * s56 + t30 * (((t29 * s16 + ddreal(6) *
     & t44) * s26 + t51) * s26 + t2 * ((s34 * t26 - s25 * (t23 * s16 + s
     &26 * t55) + t3 * (t54 - t28)) * s34 + t27) + t47 * (ddreal(28) * s
     &25 + t37)) + t31 * (t48 * s26 + t43 + t56) + t32 * ((-s25 * t55 + 
     &t49) * s34 + ddreal(6) * t28 - t52 + ddreal(24) * t47 + ddreal(10)
     & * t54) + t44 * t7 + t26 * t38 * t39 * t6 + ddreal(6) * t26 * (t26
     & * t42 + t41 * t8) + t2 * t26 * t6 * t50) * MBL1011(1) * t15
      t9 = t15 * ((t36 + (s25 * (s16 + t12) + (-s16 + t23 + s26) * s26) 
     &* s26 + s26 * (s26 - t12) * s34 + ((t11 + t23 + t25 + s56) * s56 +
     & t3 * s26 * t25 + (t21 + t11) * s26 - t10 * t9 + t17) * s56) * MBL
     &1001(0) + s25 * ((-s34 * t2 + t13) * s26 + t14) + (s25 * (t5 + t4)
     & + t6) * s26 + ((t11 + t12 + s56) * s56 + (ddreal(8) * s26 + t4) *
     & s25 - t2 * (t9 - t6)) * s56) * t1
      t7 = t18 * ((-t2 * s16 * (t26 + t6) + t28 * t3 - (t6 - s25 * (-ddr
     &eal(6) * s16 + t5)) * s26 - (-t22 * (s16 - t12 + t16) - t35) * s34
     & - ((-s34 * t3 + t2 * t34 + s56) * s56 - t2 * (t26 - t6) + ddreal(
     &6) * t24 + s34 * (s34 - s16 - t13 + t33)) * s56) * I300s16(1) + (t
     &7 - t2 * (t32 + t8) - t26 * (-t23 + t19) - (s25 * (ddreal(16) * s1
     &6 - ddreal(7) * s25) + (s25 + t19) * s26) * s26 - ((-t23 + t16) * 
     &s34 + (-s26 * t20 + t21 - t4) * t22) * s34 - ((-s34 * t20 + s25 + 
     &ddreal(9) * s26 + t19) * s56 + t3 * (s25 * s26 + s34 * t25) + ddre
     &al(9) * t6 + ddreal(16) * t24 - ddreal(7) * t26 + t23 * s34 - t29 
     &* s34) * s56 - s34 * (-t3 * (-t26 * t58 - t58 * t6) + (t22 * (-t12
     & + t34) + t35) * s34 + ((t3 * t58 + s34) * s56 - t10 * (-t59 - t60
     &) - s34 * (s34 - s16 + t5 - t16)) * s56 + t60 * t13) * MBL1110(1))
     & * t65 + t9)
      t1 = ddreal(1) / t1
      result = t7 / ddreal(4) + t65 * (((-t2 * (s25 * (s25 * t22 * t50 +
     & t62 * s25 + t63) + ((-t36 - (t46 + t26) * s34 + t62 + t64) * s56 
     &- t3 * s25 * (t47 - t62) + t26 * (t42 + t50 + t26) + t63 + t57 * (
     &-t53 + s16 * (s25 - s26))) * s56) - t26 * (-s34 * t50 + t64) - t32
     & * (-t14 - t56 + t26) - t61 * t8 + (t22)**2 * ((t5 - t4) * s25 + t
     &46) * s34 + t4 * t27) * MBL1111D2(0) + s25 * s34 * (t3 * ((t57 - t
     &59) * s25 + ((s16 + s34 - s25) * s25 - t59) * s56) - (t60 + t50) *
     & s25 - t58 * (t30 + t6)) * MBL1101(0)) * t1 * t18 + t26 * MBL1011(
     &0) * t15) / ddreal(2)

           intHLs16s25s26s34s56x1132D6eps0 = result
       end function intHLs16s25s26s34s56x1132D6eps0

       function intHLs16s25s26s34s56x1132D6eps1()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1132D6eps1
           type(dd_complex) ::  t1,t2,t3,t4

           type(dd_complex) :: result

      t1 = s26 + s56
      t2 = ddreal(2)
      t3 = ddreal(1) / ddreal(4)
      t4 = ddreal(1) / (t1)**2
      result = ddreal(1) / s16 * t4 * (t3 * (s25 * t2 - t1) + (s25)**2 *
     & MBL1011(1) / ddreal(2))

           intHLs16s25s26s34s56x1132D6eps1 = result
       end function intHLs16s25s26s34s56x1132D6eps1

       function intHLs16s25s26s34s56x1141D6eps0()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1141D6eps0
           type(dd_complex) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           type(dd_complex) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           type(dd_complex) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           type(dd_complex) ::  t131,t132,t133,t14,t15,t16,t17,t18,t19,t2,t20,t21
           type(dd_complex) ::  t22,t23,t24,t25,t26,t27,t28,t29,t3,t30,t31,t32
           type(dd_complex) ::  t33,t34,t35,t36,t37,t38,t39,t4,t40,t41,t42,t43
           type(dd_complex) ::  t44,t45,t46,t47,t48,t49,t5,t50,t51,t52,t53,t54
           type(dd_complex) ::  t55,t56,t57,t58,t59,t6,t60,t61,t62,t63,t64,t65
           type(dd_complex) ::  t66,t67,t68,t69,t7,t70,t71,t72,t73,t74,t75,t76
           type(dd_complex) ::  t77,t78,t79,t8,t80,t81,t82,t83,t84,t85,t86,t87
           type(dd_complex) ::  t88,t89,t9,t90,t91,t92,t93,t94,t95,t96,t97,t98
           type(dd_complex) ::  t99

           type(dd_complex) :: result

      t1 = ddreal(6)
      t2 = ddreal(2)
      t3 = ddreal(38)
      t4 = (s16)**2
      t5 = t2 * s25
      t6 = t1 * t4
      t7 = ddreal(17) * s25
      t8 = ddreal(4) * s25
      t9 = ddreal(12) * s16
      t10 = ddreal(5)
      t11 = t10 * s25
      t12 = t1 * s26
      t13 = (s34)**2
      t14 = (t13)**2
      t15 = s34 * t13
      t16 = s16 - s34
      t17 = s26 * t16
      t18 = ddreal(19)
      t19 = ddreal(3)
      t20 = s16 * t1
      t21 = (s26)**2
      t22 = (t21)**2
      t23 = s26 * t21
      t24 = s26 * t22
      t25 = s16 * t18
      t26 = (-t5 - t20) * s34
      t27 = t19 * (t13 + t4)
      t28 = (s56)**2
      t29 = (t28)**2
      t30 = s56 * t28
      t31 = s56 * t29
      t32 = t23 + t30
      t33 = ddreal(13) * s16
      t34 = s16 + s25
      t35 = t2 * s26
      t36 = s16 * s25
      t37 = (-s34 + t34 + s26) * s26 + (t35 + s16 + s25 - s34 + s56) * s
     &56 + t36
      t38 = ddreal(7)
      t39 = t19 * s25
      t40 = t38 * s16
      t41 = t2 * s16
      t42 = t39 + t41
      t43 = ddreal(9)
      t44 = t43 * s25
      t45 = s16 + t44
      t46 = t1 * s25
      t47 = s16 - t46
      t48 = s16 - t39
      t49 = s25 + s26
      t50 = t19 * s16
      t51 = t38 * s25
      t52 = (s25)**2
      t53 = (t52)**2
      t54 = s25 * t53
      t55 = s25 * t52
      t56 = ddreal(4) * s16
      t57 = t56 * t52
      t58 = t19 * s26
      t59 = s25 * t45
      t60 = t41 * t52
      t61 = ddreal(15) * s26
      t62 = ddreal(10) * s16
      t63 = t47 * t21
      t64 = s25 * s26
      t65 = s25 * s34
      t66 = t48 * t10
      t67 = ddreal(30)
      t68 = t67 * s26
      t69 = ddreal(12) * s25
      t70 = ddreal(25) * s25
      t71 = t10 * s16
      t72 = ddreal(18) * t21
      t73 = ddreal(8) * s16
      t74 = s16 - s26
      t75 = s16 * t52
      t76 = t66 * s26
      t77 = t43 * t55
      t78 = t19 * s34
      t79 = s25 + s34
      t80 = (t49)**2
      t81 = t49 * (((t11 - t58) * s26 + t59) * s26 + t60) * t13
      t18 = ((t28 * ((-t19 * t79 + s16) * s56 + (t41 - t69) * s25 + (t78
     & + s25 - t61 - t50) * s34 + t76) + (-t52 * t45 * t1 + ddreal(10) *
     & t48 * t21) * s26 - ddreal(12) * s25 * (-t63 + t55) + ((t65 - (s16
     & + ddreal(14) * s25 + t12) * s25 + t72) * s34 + ((-t68 - ddreal(18
     &) * s16 + t46) * s26 - t69 * (s16 - t11)) * s26 + t52 * (t71 + t70
     &)) * s34 - t73 * t55) * s56 + t52 * (t52 * (-t40 - t39) - t21 * t4
     &5 * t1) - ddreal(8) * t64 * (t42 * t52 - t63) + ((t65 * (t39 + t35
     &) - ((ddreal(28) * s26 + t44 + t50) * s25 + s26 * (t12 + t41)) * s
     &25 + ddreal(12) * t23) * s34 + t49 * (((s25 * t18 - t61 - t9) * s2
     &6 + ddreal(41) * t52) * s26 + t52 * (t62 + t44))) * s34 + t66 * t2
     &2) * s56 + t24 * t48 + t80 * (-t19 * t23 + s26 * (s25 * (t44 + t41
     &) + s26 * (t51 - t50)) + t57) * s34 - t81 + t64 * (t15 * (s26 + t3
     &9) - t55 * (t40 + t39))
      t45 = s16 * t53
      t48 = s26 + s56
      t63 = s25 + s26 - s34 + s56
      t66 = s25 + s26 + s56
      t81 = ddreal(16)
      t82 = t19 * t4
      t83 = (s16 * t81 - s25) * s25
      t84 = ddreal(11) * s16
      t85 = t58 - t5
      t86 = s16 - t5
      t87 = s16 - s25
      t88 = ddreal(23)
      t89 = t52 - t21
      t90 = t52 * s26
      t91 = ddreal(15) * s16
      t92 = s25 * t87
      t93 = t67 * t21
      t94 = -t53 + t22
      t95 = s16 * t21
      t96 = s16 + s26
      t97 = s26 * t96
      t98 = t64 * (t15 * (-s26 + t8) - t55 * (t71 + t8))
      t62 = ((t28 * (t28 + (s34 * t10 - t51 - t68 - t71) * s34 + ddreal(
     &10) * t92 + ddreal(15) * t97) + t67 * s26 * (t52 * t86 + t95) + dd
     &real(15) * t94 - ((t65 - (s16 - ddreal(13) * s25 + ddreal(12) * s2
     &6) * s25 - t93) * s34 + ((ddreal(60) * s26 + t1 * (t51 + t71)) * s
     &26 + t69 * (t50 - t8)) * s26 + t52 * (-ddreal(29) * s25 + t71)) * 
     &s34 + ddreal(60) * t92 * t21) * s56 + t13 * (-ddreal(12) * s25 * t
     &89 + (s25 * (-ddreal(26) * s25 + t41) + ddreal(20) * t21) * s26 - 
     &t75) + t2 * s34 * (t49 * (t52 * (t46 + t50) - ddreal(15) * t23 + s
     &26 * ((s25 - t62) * s26 + s25 * (s25 * t88 - t73))) + s25 * (-s26 
     &+ t5) * t13) + t23 * ((t91 + t12) * s26 + ddreal(40) * t92) - t53 
     &* (t71 + t8) - t67 * t90 * (-s26 * t86 + t52)) * s56 + t21 * t22 +
     & t80 * (((-t10 * t87 - t12) * s26 - t5 * t47) * s26 + t57) * s34 -
     & t49 * ((-t10 * t21 + ddreal(12) * t52) * s26 - t64 * t74 + t60) *
     & t13 + t98
      t68 = s16 * t54
      t73 = t87 * s26
      t88 = ddreal(11) * s25
      t98 = ddreal(21) * s16
      t99 = ddreal(14) * s16
      t100 = (s25 + t99) * s25 + ddreal(8) * t4
      t101 = (t98 + t39) * s25 + t4
      t102 = -ddreal(46) * s16
      t103 = (t102 - t44) * s25 + ddreal(15) * t4
      t104 = (s25 + t71) * s25 - ddreal(4) * t4
      t105 = t56 + t5
      t34 = s25 * t34
      t106 = t10 * t4
      t107 = -t71 - t5
      t108 = s16 + t5
      t109 = -t71 + t39
      t110 = ddreal(4) * s26
      t111 = t84 * t52
      t112 = -t103
      t113 = ddreal(8) * s25
      t114 = ddreal(33) * s16
      t115 = t114 * t53
      t116 = (ddreal(21) * s25 + t99) * s25 - t82
      t117 = -s25 + t9
      t118 = -ddreal(36) * s26
      t119 = t81 * s26
      t120 = s25 - s26
      t121 = ddreal(18) * t4
      t122 = s25 * t104
      t123 = t52 * t112
      t124 = t77 * t101
      t70 = (ddreal(82) * s16 + t70) * s25 - ddreal(11) * t4
      t125 = s25 + t41
      t126 = (-s25 + t84) * s25 + t106
      t127 = (-t51 + t56) * t125
      t128 = ddreal(10) * t126
      t129 = -ddreal(10) * t107
      t130 = ddreal(8) * s26
      t131 = ddreal(45) * t4
      t54 = ((((t98 * s26 + t131) * s26 - t122 * t67) * s26 - t123 * t1)
     & * s26 - t124) * s26 - t19 * s25 * (t15 * ((s25 - t35) * s34 - s25
     & * (t8 + t41) - (-t130 - t71 + t39) * s26) + t54) + s34 * (((s26 *
     & (-t116 * t19 + t12 * t117) + t52 * (-ddreal(54) * s16 - ddreal(18
     &) * s25)) * s25 + t129 * t23) * s34 + (t52 * t70 * t19 - t128 * t2
     &1) * s26 - ddreal(90) * s16 * t94 - ddreal(12) * s25 * (t127 * t21
     & - t53) + t121 * t55) + t45 * (-ddreal(24) * s16 - ddreal(42) * s2
     &5)
      t94 = s16 * s26
      t132 = -ddreal(63)
      t133 = s16 * t23
      t44 = ((((s16 * t132 - t44) * s25 - t82) * s25 - t110 * t112) * s2
     &5 + t14 * t2 - t93 * t104) * s25 + ((ddreal(4) * t64 * t117 + t129
     & * t21 + t39 * t4 - t99 * t52 - ddreal(21) * t55 + t65 * (-t119 - 
     &t71 + t39)) * s34 + (s25 * t70 - t130 * t127) * s25 + t21 * (-ddre
     &al(120) * t94 - t128)) * s34 + t133 * (ddreal(60) * s16 + ddreal(3
     &5) * s26)
      t70 = t4 * s25
      t70 = -t43 * t53 + ddreal(15) * s25 * (-s26 * t104 + t70) + s34 * 
     &(-t10 * t126 * s26 + ((-s25 + ddreal(10) * s26 + t9) * s25 + ddrea
     &l(25) * t94) * s34 - t70 * t81 + ddreal(14) * t55 + ddreal(20) * t
     &75 - ddreal(90) * t95) - t8 * t15 + ddreal(35) * t133 + t131 * t21
     & + t102 * t55
      t81 = -t19 * t55 + s34 * (-t107 * s34 + (t118 - t88 - t71) * s16 +
     & t52) + t69 * t4 + t121 * s26 + t98 * t21 - t91 * t52
      t93 = t1 * s34
      t99 = t64 * (-t15 * (((-t110 + t109) * s26 + t46 * t108) * s26 + t
     &111) + t45 * (t98 + t88))
      t102 = t49 * s34
      t90 = t102 * ((t19 * t64 * (s26 * ((t40 + t46) * s25 - t4) + t36 *
     & (s16 + t88)) + t23 * (t107 * s26 + s25 * (-t40 + t39)) + t6 * t55
     &) * s34 - t49 * (-t1 * t21 * (-s25 * (s25 * t105 - t4) + t95) + t2
     &3 * (-t106 + t34) + t50 * t90 * (t88 + t41) + ddreal(12) * t4 * t5
     &5))
      t95 = t21 * ((-t101 * t55 + t21 * (s26 * t4 - t122)) * s26 - t100 
     &* t53)
      t101 = t94 + t52
      t104 = t52 + t4
      t87 = (t87)**2
      t106 = t87 * s26
      t112 = s25 * t104
      t79 = t2 * s16 * t79
      t126 = t104 * t52
      t127 = t87 * t21
      t128 = s25 * (-t106 + t75)
      t129 = t49 * s26
      t130 = ddreal(1) / s16
      t63 = ddreal(1) / (t63)**2
      t48 = ddreal(1) / (t48)**2
      t37 = ddreal(1) / (t37)**2
      t22 = (t19 * t95 + s16 * t23 * t22 + s16 * t30 * t29 + (-t1 * t64 
     &* (t100 * t55 + t120 * t14) + t21 * ((((t40 * s26 + t121) * s26 - 
     &ddreal(15) * t122) * s26 - ddreal(4) * t123) * s26 - t124) + ((t65
     & * (((t19 * t109 - t119) * s26 + t69 * t108) * s26 + t111) + t8 * 
     &t117 * t23 - t115 + t118 * t55 * (s25 + t50) - t39 * t116 * t21 - 
     &t107 * t10 * t22 - t77 * t4) * s34 + t49 * (t19 * t64 * (s25 * ((d
     &dreal(49) * s16 + t113) * s25 + t2 * t4) + s26 * ((t114 + t7) * s2
     &5 - ddreal(13) * t4)) + t23 * ((t11 - t25) * s25 - ddreal(25) * t4
     &) + t115 - ddreal(36) * s16 * t22 + t67 * t4 * t55)) * s34 - t68 *
     & (t98 + t88)) * s56 + t28 * t54 + t29 * t70 + t30 * t44 + t31 * t8
     &1 - s25 * (t39 - t35) * t21 * t14 + s16 * (s26 * t38 + t50 - t93) 
     &* t28 * t29 + t52 * t103 * t22 - t99 - t90 - t6 * t52 * t53) * MBL
     &1011(1)
      t14 = (t2 * t87 * t49 * t80 * s34 + (-t19 * t14 * (s34 - s25 - t11
     &0 - t41) - t13 * (((t51 + t50 + t35) * s16 - t43 * t89 + ddreal(12
     &) * t64) * s34 + t49 * ((-t58 - t33) * s25 + ddreal(15) * t101 + d
     &dreal(10) * t4)) + t93 * t87 * t80) * s56 + t14 * (-t85 * s34 + (s
     &16 + t58 - t113) * s25 + t1 * t97) + t28 * (t1 * s34 * (t49 * t87 
     &+ t15) - t13 * ((s26 * t43 + s16 + t46) * s34 + s16 * (t61 + t71) 
     &+ (s16 + t46 - t58) * s25)) - s34 * ((-s25 + t71 + t78) * s34 - t1
     &04 * t2 + t56 * s25) * t30 + t49 * ((-t20 + t69) * s25 - t19 * (t1
     &29 + t4) - t94) * t15 - t80 * (t10 * s16 * t96 + (-s16 * t43 - s26
     & + t113) * s25) * t13) * MBL1110(1) * t130
      t10 = t48 * ((t19 * (s16 * t24 + t31 * (t2 * (s26 - s34) + s16)) -
     & t2 * (-t30 * ((t91 * s26 + ddreal(20) * t92) * s26 - t2 * s34 * (
     &-(s26 * t10 + s25) * s34 + s25 * (t50 - t8) + (t61 + t51 + t71) * 
     &s26) + ddreal(10) * t23 - ddreal(10) * t55 + t71 * t52) + t68) + t
     &62 + ddreal(10) * s25 * t23 * (s25 * t86 + t73) - ddreal(15) * t53
     & * t21) * MBL1001(0) - t2 * (s25 * (t23 * (-s26 * t47 + t59) + t45
     &) + t30 * ((-t8 * t47 - t76) * s26 + ((s25 - t12) * s34 + t2 * s25
     & * t74 + (t61 + t20) * s26 - ddreal(10) * t52) * s34 + t75 + t77))
     & + t18 - ddreal(4) * t55 * t42 * t21 + t22 * t130) + t14
      t14 = ddreal(1) / t66
      t4 = t14 * (t130 * ((t65 * (-t2 * (-(-t60 + (s26 * s34 - t41 * s26
     & - t34) * s34 + t106 + t112) * s56 + t128 + t102 * t101) + t13 * (
     &t52 + t21) + t28 * (-t79 + t52 + t13 + t4) + t126 + t127) * MBL110
     &1(0) + (-s34 * t120 - (s16 - s25 - s34) * s56 - t73 - t92) * (-t2 
     &* t128 + ((t129 + t52) * s34 - t49 * (s26 * t125 + t2 * t52 + t36)
     &) * s34 + ((-t79 + (s34 - s25) * s34 + t4 + t52) * s56 + t2 * (t11
     &2 + t106) - t57 - s34 * (s26 * t105 - (s25 + t35) * s34 + t19 * t3
     &4)) * s56 + t126 + t127) * MBL1111D2(0)) * t14 * t37 - t15 * MBL11
     &10(0) * t63 - t55 * MBL1011(0) * t48) + I300s16(0))
      result = t37 * ((s25 * ((-t11 + t33) * s25 + t6) + ddreal(14) * t3
     &2 + ((ddreal(25) * s16 + t7) * s26 + (s16 * t3 - t5) * s25 + t6) *
     & s26 + ((t12 - t11) * s34 - s25 * (s16 - ddreal(10) * s25) - (ddre
     &al(25) * s26 + t8 + t9) * s26) * s34 + (ddreal(42) * s26 + ddreal(
     &25) * t16 + t7) * t28 + t2 * (s25 * (-s25 + t25) + ddreal(25) * t1
     &7 + ddreal(21) * t21 + t26 + t27 + t7 * s26) * s56) * t130 + (s25 
     &* ((t71 - t5) * s25 + t82) + t1 * t32 + ((t84 + t51) * s26 + t82 +
     & t83) * s26 + (t85 * s34 - s25 * (s16 - t8) - (ddreal(11) * s26 + 
     &t20 + t5) * s26) * s34 + ((ddreal(18) * s26 + ddreal(11) * t16 + t
     &51) * s56 + ddreal(22) * t17 + t26 + t27 + ddreal(14) * t64 + t72 
     &+ t83) * s56) * I300s16(1) + t10 * t63) / ddreal(12) + t4 / ddreal
     &(6)

           intHLs16s25s26s34s56x1141D6eps0 = result
       end function intHLs16s25s26s34s56x1141D6eps0

       function intHLs16s25s26s34s56x1141D6eps1()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1141D6eps1
           type(dd_complex) ::  t1,t10,t11,t12,t13,t2,t3,t4,t5,t6,t7,t8
           type(dd_complex) ::  t9

           type(dd_complex) :: result

      t1 = -ddreal(11)
      t2 = (s25)**2
      t3 = ddreal(2) * s25
      t4 = ddreal(3)
      t5 = -ddreal(14)
      t6 = ddreal(18) * s26
      t7 = (s26)**2
      t8 = ddreal(6)
      t9 = s26 + s56
      t10 = s25 + s26 - s34 + s56
      t11 = s25 + s26 + s56
      t11 = ddreal(1) / t11
      t12 = ddreal(1) / s16
      t10 = ddreal(1) / (t10)**2
      t9 = ddreal(1) / (t9)**2
      t13 = ddreal(1) / ddreal(12)
      result = t11 * (-((s34)**3 * MBL1110(1) * t10 + s25 * t2 * MBL1011
     &(1) * t9) * t12 + I300s16(1)) / ddreal(6) + t13 * ((-ddreal(2) * t
     &2 + ddreal(7) * t7) * s25 + t8 * ((s56)**3 + s26 * t7) - s26 * t2 
     &+ ((s26 * t4 - t3) * s34 + (s26 * t1 - t3) * s26 + ddreal(4) * t2)
     & * s34 + ((s34 * t1 + ddreal(7) * s25 + t6) * s56 - (s25 * t5 - t6
     &) * s26 - (-s34 * t4 + ddreal(22) * s26 + t3) * s34 - t2) * s56) *
     & t12 * t9 * t10

           intHLs16s25s26s34s56x1141D6eps1 = result
       end function intHLs16s25s26s34s56x1141D6eps1

       function intHLs16s25s26s34s56x1211D2eps0()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1211D2eps0
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s25 + s26 + s56
      t1 = ddreal(1) / t1
      result = -ddreal(1) / s16 * ((s25 * MBL1101(0) + (s26 + s56) * MBL
     &1111D2(0)) * t1 - MBL1110(0))

           intHLs16s25s26s34s56x1211D2eps0 = result
       end function intHLs16s25s26s34s56x1211D2eps0

       function intHLs16s25s26s34s56x1211D2eps1()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1211D2eps1

           type(dd_complex) :: result

      result = MBL1110(1) / s16

           intHLs16s25s26s34s56x1211D2eps1 = result
       end function intHLs16s25s26s34s56x1211D2eps1

       function intHLs16s25s26s34s56x1211D4eps0()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1211D4eps0
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           type(dd_complex) ::  t20,t21,t22,t23,t24,t3,t4,t5,t6,t7,t8,t9

           type(dd_complex) :: result

      t1 = s16 - s25
      t2 = ddreal(2) * s16
      t3 = s25 + t2
      t4 = s16 * s25
      t5 = (s25)**2
      t6 = (s16)**2
      t7 = (s34)**2
      t8 = t6 * s25
      t9 = s26 * t7
      t10 = (s26 * t1 - t5) * s16
      t11 = (s26 * t3 + t4) * s34
      t12 = t5 + t7 + t6
      t13 = s16 + s25
      t14 = -ddreal(2) * s34 * t13 + t12 - ddreal(2) * t4
      t15 = s16 + s25 - s34
      t16 = (-s34 + t13 + s26) * s26 + (ddreal(2) * s26 + t15 + s56) * s
     &56 + t4
      t17 = s26 + s56
      t18 = s26 * s34
      t19 = -ddreal(3) * s16
      t20 = ddreal(2) * s25
      t21 = s25 - t2
      t22 = -ddreal(4)
      t23 = ddreal(1) / t16
      t16 = ddreal(1) / t16
      t14 = ddreal(1) / t14
      t24 = ddreal(1) / ddreal(2)
      result = t24 * (t17 * (s25 + s26 + s56) * MBL1111D2(0) * t16 + (((
     &((t19 - t20) * s25 + t6) * s26 + (s26 * t21 + t18 - t4) * s34 + ((
     &s26 * t22 + t19 - t20) * s25 + (t21 + s34) * s34 + t6) * s56 + t8 
     &- ddreal(2) * s25 * ((s26)**2 + (s56)**2) + t19 * t5) * MBL1110(1)
     & - s25 * (t2 * s25 - (s34 - t13) * s26 + s56 * t15) * MBL1101(0)) 
     &* s34 + (-s16 * I300s16(1) + MBL1001(0)) * ((-s25 * s34 - t2 * s34
     & - t4 + t6 + t7) * s56 + t10 - t11 + t8 + t9) + t17 * (-t2 * t5 + 
     &(t5 + t6) * s25 + (t1)**2 * s26 + (-s25 * t13 - t2 * s26 + t18) * 
     &s34 + (t12 - ddreal(2) * s16 * (s25 + s34)) * s56) * MBL1011(1)) *
     & t14 * t23) - ((s16 * t1 + (-t3 + s34) * s34) * s56 + t8 + t9 + t1
     &0 - t11) * t14 * t23

           intHLs16s25s26s34s56x1211D4eps0 = result
       end function intHLs16s25s26s34s56x1211D4eps0

       function intHLs16s25s26s34s56x1211D4eps1()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1211D4eps1

           type(dd_complex) :: result

      result = ddreal(0)

           intHLs16s25s26s34s56x1211D4eps1 = result
       end function intHLs16s25s26s34s56x1211D4eps1

       function intHLs16s25s26s34s56x1212D4eps0()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1212D4eps0
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t2,t3,t4,t5,t6,t7
           type(dd_complex) ::  t8,t9

           type(dd_complex) :: result

      t1 = s16 - s25
      t2 = s25 * t1 + s26 * t1 + (s25 - s26) * s34 + (s16 - s25 - s34) *
     & s56
      t3 = s16 + s25
      t4 = s16 * s25
      t5 = (s16)**2
      t6 = (s25)**2
      t7 = (s34)**2 + t5 + t6
      t8 = ddreal(2)
      t9 = -t8 * (s34 * t3 + t4) + t7
      t10 = s16 + s25 - s34
      t11 = (-s34 + t3 + s26) * s26 + (s26 * t8 + s56 + t10) * s56 + t4
      t12 = s26 + s56
      t13 = t8 * t6
      t14 = s25 * s34
      t4 = t4 * t8
      t9 = ddreal(1) / t9
      t11 = ddreal(1) / t11
      result = t8 * t2 * t9 * t11 - t11 * (-(s25 + s26 + s56) * MBL1111D
     &2(0) + t9 * (t2 * MBL1001(0) + (-s34 * (s25 + s26 - s34 + s56) * (
     &t4 - (s34 - t3) * s26 + s56 * t10) * MBL1110(1) + t12 * (-t13 * s1
     &6 + s25 * t7 + (t1)**2 * s26 - (s26 * t3 + t13) * s34 + (-s16 * s3
     &4 - t14 - t4 + t5 + t6) * s56) * MBL1011(1)) / s16 - s16 * t2 * I3
     &00s16(1) - t14 * (t12 * t8 + t10) * MBL1101(0)))

           intHLs16s25s26s34s56x1212D4eps0 = result
       end function intHLs16s25s26s34s56x1212D4eps0

       function intHLs16s25s26s34s56x1212D4eps1()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1212D4eps1

           type(dd_complex) :: result

      result = ddreal(0)

           intHLs16s25s26s34s56x1212D4eps1 = result
       end function intHLs16s25s26s34s56x1212D4eps1

       function intHLs16s25s26s34s56x1213D6eps0()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1213D6eps0
           type(dd_complex) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           type(dd_complex) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           type(dd_complex) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           type(dd_complex) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           type(dd_complex) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           type(dd_complex) ::  t153,t154,t16,t17,t18,t19,t2,t20,t21,t22,t23,t24
           type(dd_complex) ::  t25,t26,t27,t28,t29,t3,t30,t31,t32,t33,t34,t35
           type(dd_complex) ::  t36,t37,t38,t39,t4,t40,t41,t42,t43,t44,t45,t46
           type(dd_complex) ::  t47,t48,t49,t5,t50,t51,t52,t53,t54,t55,t56,t57
           type(dd_complex) ::  t58,t59,t6,t60,t61,t62,t63,t64,t65,t66,t67,t68
           type(dd_complex) ::  t69,t7,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79
           type(dd_complex) ::  t8,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t9
           type(dd_complex) ::  t90,t91,t92,t93,t94,t95,t96,t97,t98,t99

           type(dd_complex) :: result

      t1 = s16 - s25
      t2 = ddreal(11)
      t3 = (s16)**2
      t4 = (t3)**2
      t5 = s16 * t4
      t6 = s16 * t3
      t7 = ddreal(6) * t3
      t8 = ddreal(3)
      t9 = ddreal(2)
      t10 = ddreal(7) * s16
      t11 = t9 * s25
      t12 = t8 * t3
      t13 = (-t11 - t10) * s25 + t12
      t14 = s16 + s25
      t15 = s25 + s26
      t16 = t8 * s16
      t17 = ddreal(5)
      t18 = t17 * s25
      t19 = t9 * t3
      t20 = ddreal(12) * s16
      t21 = t8 * s25
      t22 = t17 * t3
      t23 = (t21 + t20) * s25 + t22
      t24 = (s25)**2
      t25 = (t24)**2
      t26 = s25 * t24
      t27 = t26 * t25
      t28 = s25 * t25
      t29 = ddreal(20) * s16
      t30 = ddreal(10) * s25
      t31 = ddreal(63) * s16
      t32 = ddreal(30) * t3
      t33 = t17 * t6
      t34 = ddreal(6) * s25
      t35 = (ddreal(23) * s16 + t34) * s25 + t3
      t36 = ddreal(8) * s16
      t37 = ddreal(7) * t3
      t38 = (t21 + t36) * s25 - t37
      t39 = (s26)**2
      t40 = s26 * t39
      t41 = ddreal(30) * s25
      t42 = s16 * s25
      t43 = ddreal(18)
      t44 = s16 * t43
      t45 = ddreal(8) * s25
      t46 = (-t44 - t45) * s25 + t7
      t47 = t16 - s25
      t48 = s25 * ((ddreal(47) * s16 + t30) * s25 + ddreal(17) * t3)
      t49 = s16 * t24
      t50 = ddreal(13) * s16
      t51 = -t16 + t11
      t52 = ddreal(4)
      t53 = ddreal(17) * s16
      t54 = t52 * s25
      t55 = ddreal(12) * t3
      t56 = t8 * t6
      t57 = ddreal(42)
      t58 = ddreal(6) * t6
      t59 = t8 * s26
      t60 = t59 * t23
      t61 = t8 * t38
      t62 = t52 * s26
      t63 = ddreal(6) * s16
      t64 = t63 * s26
      t65 = s25 * s34
      t66 = (t1)**2
      t67 = (t66)**2
      t68 = t1 * t67
      t69 = t1 * t66
      t70 = t59 * t14
      t71 = t4 * t24
      t72 = t69 * s26
      t73 = -ddreal(34)
      t74 = ddreal(24) * t3
      t75 = ddreal(12) * s25
      t76 = ddreal(9) * s26
      t77 = s16 * t25
      t78 = -t3 + t24
      t79 = t1 * t15
      t80 = t14 * s26
      t73 = -(-t9 * t42 * t78 - ((-s34 * t47 - t38) * s34 - t1 * t23) * 
     &s34 + t25 - t4) * s56 - s25 * (((t24 * t9 - t3 * t43) * s25 + ddre
     &al(28) * t6) * s25 - ddreal(16) * t4) - t8 * (-t72 * t14 + t5) + (
     &((-s34 * t51 - s16 * (t63 + t76) - (-t44 - t45 - t59) * s25) * s34
     & - s25 * ((ddreal(46) * s16 + t75) * s25 + t19) - s26 * t61) * s34
     & + t1 * (((s16 * t73 - t45) * s25 - t74) * s25 + t58 - t60)) * s34
     & - t77
      t60 = t73 * s56 - t36 * t28 - ((t54 * t35 * s26 + t33 * s25 + (((t
     &50 + t62 + t18) * s25 - t64 - t65) * s34 - t48 - s26 * (-t8 * t47 
     &* s26 - t9 * t46)) * s34 + t32 * t24 + t31 * t26 + t39 * t61 + ddr
     &eal(10) * t25) * s34 - t1 * (s25 * (((-s16 * t57 - t18) * s25 - dd
     &real(19) * t3) * s25 + t58) + (t52 * (((-t54 - t53) * s25 - t55) *
     & s25 + t56) - t60) * s26)) * s34 - t24 * t25 + ddreal(29) * t71 + 
     &ddreal(36) * t3 * t25 - ddreal(50) * t6 * t26 - t34 * t5 - t72 * (
     &t13 * t9 - t70)
      t61 = (s34)**2
      t73 = (t61)**2
      t81 = s34 * t73
      t82 = s34 * t61
      t83 = t61 + t3 + t24
      t84 = s34 * t14
      t85 = -t9 * (t84 + t42) + t83
      t86 = s16 + s25 - s34
      t87 = t9 * s26
      t88 = (-s34 + t14 + s26) * s26 + (t87 + t86 + s56) * s56 + t42
      t89 = t9 * s16
      t90 = -t54 + s16
      t91 = ddreal(9) * s16
      t92 = t52 * s16
      t93 = t34 + t92
      t94 = t42 * t14
      t95 = ddreal(8) * t26 + (t93 * s26 + (t91 + t75) * s25 - t3) * s26
     & + t94
      t96 = ddreal(15) * s25
      t97 = t11 + s16
      t98 = ddreal(22) * s25
      t99 = (t98 - s16) * s25 - t12
      t100 = t17 * s16
      t101 = -t100 + t21
      t102 = t24 * ((t75 + t10) * s25 + t12) + ((s26 * t101 + t99) * s26
     & - s25 * (-t96 + s16) * t97) * s26
      t103 = t34 - t100
      t104 = t103 * t39 + t24 * (t45 + t10) - t9 * s26 * (-s25 * (t30 + 
     &t16) + t39)
      t105 = t9 * (-t39 + t24) + t18 * s26
      t99 = t9 * t99
      t106 = t8 * t101
      t107 = t3 * s25
      t108 = (-t50 - t41) * t24
      t109 = t62 - t18
      t110 = s34 * (t24 * (-t30 - t100) + s26 * (-ddreal(3) / ddreal(2) 
     &* t93 * s26 + (-s25 * t43 - t50) * s25 + t3))
      t111 = s26 * t90 * t66
      t103 = (s26 * t103 - t8 * (-t42 + t39) + ddreal(10) * t24) * t82
      t112 = t4 * s25
      t113 = ddreal(21) * t3
      t114 = (-t18 + t53) * s25
      t115 = t2 * t6
      t116 = t69 * t39
      t117 = t116 * t8
      t118 = ddreal(12) * s26
      t119 = ddreal(13) * s25
      t120 = s25 * (((-t50 + t54) * s25 + ddreal(15) * t3) * s25 - ddrea
     &l(7) * t6) - t73 * t9 + (((-ddreal(6) * s25 + ddreal(6) * s26 + t1
     &00) * s34 - s16 * (ddreal(15) * s26 + t16) - (-t98 - t76 + s16) * 
     &s25) * s34 - t1 * (-t43 * s25 * t15 + (-t119 - t118 + s16) * s16))
     & * s34 + t4 - t59 * t69
      t93 = t9 * s34 * (t1 * t93 / ddreal(2) + t61)
      t121 = t8 * t42 * t1
      t122 = (s56)**2
      t123 = s56 * t122
      t124 = s25 * (-t18 + t89)
      t125 = ddreal(9) * t3
      t126 = (-t11 + t10) * s25
      t127 = s16 * t26
      t128 = t26 * ((-t126 + t125) * s25 - t33) + (((t61 * t101 + t121 +
     & t26 - t6 + t93) * s56 + t120) * s56 + t24 * ((-t114 + t113) * s25
     & - t115) + t9 * (t1 * (t111 - t110) - t103 + t112) - t117 + t61 * 
     &(-t61 * t109 - t108 + (t106 * s26 + t99) * s26 - t107)) * s56 + t1
     &02 * t61 - t104 * t82 + t105 * t73 - t40 * t69 + t71 + t116 * t90 
     &+ t79 * t95 * s34 + t124 * t72
      t129 = -t21 + s16
      t130 = t14 * t39
      t131 = t3 * t52
      t132 = t3 + t24
      t133 = ddreal(6) * t42
      t134 = -t132 * t17 + t133
      t135 = -t17 * s16 * t14 + ddreal(12) * t24
      t136 = (t96 - t63) * s25
      t137 = t136 + t22
      t138 = t2 * s25
      t139 = -t138 + t89
      t140 = t9 * t97
      t141 = t87 + s25
      t142 = ddreal(9) * t6
      t143 = ((ddreal(51) * s16 - ddreal(63) * s25) * s25 + ddreal(33) *
     & t3) * s25 - t142
      t144 = -t134 * t52
      t145 = t6 * t43
      t142 = s25 * (((ddreal(29) * s16 - ddreal(26) * s25) * s25 + t74) 
     &* s25 + t142)
      t146 = t26 * t1
      t147 = t137 * t8
      t148 = t129 * t67
      t149 = ddreal(7) * s25
      t150 = ddreal(22) * t3
      t151 = t52 * t67
      t119 = (-t92 - t119) * s25
      t152 = t6 - t26
      t153 = t146 * t3
      t154 = (t15)**2
      t53 = ((t134 * t61 + t52 * (t42 * t132 + t84 * t66) - t25 - t4 + t
     &9 * t14 * t82 - t7 * t24) * s56 + s25 * ((((t50 - t21) * s25 - t15
     &0) * s25 + t145) * s25 - ddreal(7) * t4) + (((t22 + ddreal(8) * t8
     &0 + t136) * s34 + ((-ddreal(21) * s25 + t53) * s25 + t2 * t3) * s2
     &5 - t56 - t144 * s26) * s34 - t66 * (t3 - ddreal(16) * t80 + t119)
     &) * s34 + t5 - t140 * t73 - t151 * s26) * s56 - t8 * ((-t5 + t28) 
     &* s25 - t148 * s26) - (((-((-t138 + t89 - t118) * s25 - t64) * s34
     & - (t147 + ddreal(12) * t80) * s26 - ddreal(24) * t26 + ddreal(10)
     & * t94) * s34 - (ddreal(6) * t134 * s26 + t143) * s26 - t142) * s3
     &4 + t66 * (s25 * ((-ddreal(14) * s25 + t36) * s25 + t7) + (t8 * (t
     &119 + t3) - ddreal(24) * t80) * s26)) * s34 - ddreal(15) * t49 * t
     &152 + ddreal(30) * t153 - ddreal(6) * t67 * t39 + t11 * t81
      t20 = t53 * s56 + t8 * (t148 * t39 + t24 * t5) - ((((-(s25 * (-t18
     & + s16) + t87 * t139) * s25 + ddreal(6) * t97 * t39 - t65 * (t62 +
     & s25)) * s34 - ((t147 + ddreal(8) * t80) * s26 + t54 * t135) * s26
     & + ddreal(10) * t146 + ddreal(14) * t3 * t24) * s34 - t142 * t87 -
     & t143 * t39 + t144 * t40 - t145 * t24 - t55 * t26 + ddreal(10) * t
     &28 - ddreal(24) * t77) * s34 + t66 * t15 * (s25 * ((-t18 + t20) * 
     &s25 + t125) + ((-ddreal(23) * s25 + t92) * s25 + t12 - ddreal(16) 
     &* t80) * s26)) * s34 - t27 - t151 * t40 + t34 * t68 * s26 + t127 *
     & (((t149 - t44) * s25 + t150) * s25 - ddreal(13) * t6)
      t5 = (t39 * t68 * t8 + t58 * t25) * s25 - t52 * t3 * t25 * t132 + 
     &s16 * t27 + t20 * s56 + t26 * t5 + t73 * (s25 * ((s25 * (-t18 + s1
     &6) + t139 * s26) * s26 + t49) - t140 * t40) + t82 * (t137 * t40 - 
     &t52 * t127 * t14 + t9 * s26 * ((s25 * t135 + t130) * s26 + t24 * (
     &-t17 * s25 * t1 - t37))) + t15 * ((t131 * t24 + t87 * (((t91 - t18
     &) * s25 + t131) * s25 + t58)) * s25 + t39 * (t134 * s26 + (t2 * t2
     &4 - t12) * s16 - ddreal(16) * s25 * t78) + ddreal(6) * t49 * t132)
     & * t61 + s25 * s26 * t141 * t81 - t66 * t154 * (-t52 * (t130 - t94
     &) + t1 * (t18 + s16) * s26) * s34 + t67 * s26 * (t24 * t47 + t39 *
     & (t129 - s26))
      t20 = (t54 + t100) * s25
      t25 = t17 * s26
      t27 = (t149 + t36) * s25
      t28 = s25 + s26 - s34 + s56
      t37 = ddreal(9) * s25
      t40 = ddreal(10) * t42
      t53 = t40 + t132
      t43 = t42 * t43 + t78 * t8
      t55 = t42 * (t37 + t100)
      t56 = t6 - t26
      t11 = (((-s34 * t39 + t39 * (s26 + t21 + s16) + t42 * (ddreal(7) *
     & s26 + t11)) * s34 - t39 * t43 - t9 * s26 * (t130 + t55) - t49 * (
     &t34 + t36)) * s34 + t123 * (-t84 * t9 + t40 + t83) + t15 * (s26 * 
     &(t53 * s26 + s16 * ((t29 + t37) * s25 - t22)) + t42 * ((t34 + t89)
     & * s25 + t131))) * s34 + s56 * ((t9 * s16 * t152 - (t61 * (-t15 * 
     &t8 - s16 + s34) - ((ddreal(19) * s16 + s25) * s25 + t113) * s25 + 
     &t33 - t59 * t53) * s34 - ddreal(6) * t107 * t1 + t8 * (-t80 * t9 -
     & t133 - t78) * t61) * s56 + t52 * s16 * (s25 * t56 + t72) - t9 * t
     &61 * (s26 * t61 + (t70 + t43) * s26 + t55) + s34 * (t61 * ((t34 + 
     &t89 + t59) * s26 + ddreal(7) * t42) - s25 * t6 + (((ddreal(38) * s
     &16 + t11) * s25 + t3 * t57) * s25 - ddreal(10) * t6) * s26 + t8 * 
     &t53 * t39 + t49 * (ddreal(22) * s16 + t96)) - ddreal(12) * t3 * t2
     &4 * t1)
      t34 = s25 + s26 + s56
      t36 = ddreal(1) / (t88)**2
      t37 = ddreal(1) / (t85)**2
      t2 = -t5 * MBL1011(1) - ((((-t49 * t17 - s26 * (s25 * (t50 + t18) 
     &+ s26 * t51) + t65 * s26) * s34 - t49 * (-ddreal(16) * s16 - ddrea
     &l(20) * s25) + ((-s26 * t47 - t46) * s26 + t48) * s26) * s34 - s25
     & * (t9 * t35 * t39 + t42 * ((ddreal(14) * s16 + t41) * s25 + ddrea
     &l(8) * t3)) - s26 * (s25 * (((t31 + t30) * s25 + t32) * s25 + t33)
     & + t38 * t39)) * s34 + t79 * (-t29 * t26 + s26 * (-t23 * s26 + (t1
     &6 + s25) * ((-t18 - t10) * s25 + t19)))) * s34 - t60 * s56 - t72 *
     & (-s25 * ((-s16 * t2 - s25) * s25 + t7) + (t80 - t13) * s26) + t49
     & * ((((t18 - t44) * s25 + t74) * s25 - ddreal(14) * t6) * s25 + t4
     & * t8) + s34 * t28 * (t9 * s16 * (t24 * t56 + t116) + t11 - ddreal
     &(6) * t153 + t72 * t42 * t52) * MBL1110(1)
      t3 = t36 * ((t65 * ((((s34 * t141 - t52 * s25 * t14 - (t25 + t63 +
     & t45) * s26) * s34 + t17 * t94 + t39 * (t138 + t91) + t8 * s26 * (
     &t39 + t20 + t19) + ddreal(6) * t26) * s34 - t15 * ((t70 + t20 + t1
     &2) * s26 + t26 * t52 + t9 * s16 * ((s16 - s25) * s25 + t3))) * s34
     & + ((-t17 * t82 + t121 + ((ddreal(9) * s16 + ddreal(9) * s26 + t13
     &8) * s34 - t12 - t27 - ddreal(9) * t80) * s34 + t26 - t6) * s56 + 
     &t9 * (-s25 * t152 + t82 * (-t54 - t25 - t16 + s34) - t72) - s34 * 
     &(-((ddreal(15) * s16 + t75) * s25 + (t98 + t44 + t76) * s26 + t7) 
     &* s34 + s25 * ((t16 + t45) * s25 + t22) + t9 * ((t27 + t12) * s26 
     &+ t6) + ddreal(9) * t130) + ddreal(6) * t49 * t1) * s56 - t154 * t
     &69) - t21 * t86 * t61 * t123) * MBL1101(0) * t37 + (t34)**2 * t28 
     &* MBL1111D2(0))
      result = -t36 * t37 * (-(-s16 * (t61 * ((s34 * t105 - t104) * s34 
     &+ t102) + t71) + (s16 * (t24 * ((t114 - t113) * s25 + t115) + t117
     &) + t9 * s16 * (t1 * (-t111 + t110) + t103 - t112) + t61 * (s16 * 
     &(t108 + (-t106 * s26 - t99) * s26 + t107) + t61 * s16 * t109)) * s
     &56 + s16 * (-t61 * t101 - t121 - t26 + t6 - t93) * t123 - s16 * t1
     &20 * t122 - t72 * s16 * ((t90 - s26) * s26 + t124) - t79 * s16 * t
     &95 * s34 + t127 * ((t126 - t125) * s25 + t33)) * I300s16(1) - t128
     & * MBL1001(0) + t2 / s16) / ddreal(4) + t3 / ddreal(2)

           intHLs16s25s26s34s56x1213D6eps0 = result
       end function intHLs16s25s26s34s56x1213D6eps0

       function intHLs16s25s26s34s56x1213D6eps1()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1213D6eps1

           type(dd_complex) :: result

      result = ddreal(0)

           intHLs16s25s26s34s56x1213D6eps1 = result
       end function intHLs16s25s26s34s56x1213D6eps1

       function intHLs16s25s26s34s56x1221D4eps0()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1221D4eps0
           type(dd_complex) ::  t1,t10,t11,t2,t3,t4,t5,t6,t7,t8,t9

           type(dd_complex) :: result

      t1 = s26 + s56
      t2 = ddreal(2)
      t3 = s16 * s25
      t4 = (-s34 + s16 + s25 + s26) * s26 + (s26 * t2 + s16 + s25 - s34 
     &+ s56) * s56 + t3
      t5 = s25 + s26 + s56
      t6 = s25 + s26 - s34 + s56
      t7 = t2 * s16
      t8 = s16 - s34
      t4 = ddreal(1) / t4
      t9 = ddreal(1) / t5
      t6 = ddreal(1) / t6
      t10 = ddreal(1) / s16
      t11 = t4 * t10
      result = -t11 * t2 * t1 - t6 * (t10 * (s34 * MBL1110(0) + (s34 * (
     &t7 * s25 + (-s34 + s25 + t7 + s26) * s26 + (t2 * (s26 + s16) + s25
     & - s34 + s56) * s56) * MBL1110(1) + t1 * ((s16 + s34 - s25) * s25 
     &+ t1 * (s16 - s25)) * MBL1011(1)) * t4) + t5 * MBL1001(0) * t4) - 
     &t1 * I300s16(1) * t4 + t11 * ((-t2 * ((s26 * t8 + t3) * s56 + t3 *
     & s26) - s16 * (s25)**2 - t8 * (s26)**2 - t8 * (s56)**2) * MBL1111D
     &2(0) + s25 * s34 * t1 * MBL1101(0)) * t9

           intHLs16s25s26s34s56x1221D4eps0 = result
       end function intHLs16s25s26s34s56x1221D4eps0

       function intHLs16s25s26s34s56x1221D4eps1()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1221D4eps1
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s25 + s26 - s34 + s56
      t1 = ddreal(1) / t1
      result = -ddreal(1) / s16 * t1 * (s34 * MBL1110(1) + ddreal(1))

           intHLs16s25s26s34s56x1221D4eps1 = result
       end function intHLs16s25s26s34s56x1221D4eps1

       function intHLs16s25s26s34s56x1222D6eps0()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1222D6eps0
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           type(dd_complex) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           type(dd_complex) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           type(dd_complex) ::  t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52
           type(dd_complex) ::  t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63
           type(dd_complex) ::  t64,t65,t66,t67,t7,t8,t9

           type(dd_complex) :: result

      t1 = s16 - s25
      t2 = ddreal(5)
      t3 = ddreal(2)
      t4 = t3 * s25
      t5 = t2 * s16
      t6 = t5 + t4
      t7 = s25 + s26
      t8 = ddreal(6)
      t9 = s16 * t8
      t10 = ddreal(3) * s25
      t11 = t10 + t9
      t12 = (s16)**2
      t13 = s16 * t12
      t14 = ddreal(11) * t12
      t15 = s16 + s25
      t16 = ddreal(10) * s25
      t17 = t15 * s26
      t18 = s16 * s25
      t19 = ddreal(8) * s16
      t20 = ddreal(4)
      t21 = t20 * s25
      t22 = ddreal(7) * s16
      t23 = (s25)**2
      t24 = (t23)**2
      t25 = s25 * t24
      t26 = s25 * t23
      t27 = t5 * t23
      t28 = t2 * s25
      t29 = (t28 + t19) * s25
      t30 = (s26)**2
      t31 = t3 * s26
      t32 = (s25 + t31) * s34
      t33 = (t1)**2
      t34 = ddreal(3) * s26
      t35 = t33 * s26
      t36 = s16 * t23
      t37 = t3 * t23
      t38 = t2 * t13
      t39 = ddreal(3) * t35
      t40 = (s34)**2
      t41 = s34 * t40
      t42 = -t3 * (s34 * t15 + t18) + t12 + t23 + t40
      t43 = s26 * t7
      t44 = t12 + t23
      t45 = (-s34 + t15 + s26) * s26 + (t31 + s16 + s25 - s34 + s56) * s
     &56 + t18
      t46 = t3 * s16
      t47 = s25 + t46
      t48 = t47 * s26
      t49 = s25 * t15
      t50 = t20 * s16
      t51 = t36 * t20
      t52 = t3 * (s25 * t44 + t35) - s34 * (-t32 + (t50 + t4) * s26 + dd
     &real(3) * t49) - t51
      t53 = (-t3 * s16 * (s25 + s34) + (s34 - s25) * s34 + t12 + t23) * 
     &s56
      t54 = t33 * t30
      t55 = t44 * t23
      t56 = (-(t43 + t23) * s34 + t7 * (t37 + t48 + t18)) * s34
      t57 = t3 * s25 * (-t36 + t35)
      t58 = t5 - t4
      t59 = -s25 + t46
      t60 = t20 * t12
      t61 = -s25 * t1 - t60
      t62 = ddreal(3) * s16
      t63 = t23 * t1
      t64 = s16 * t15
      t65 = ddreal(9) * s16
      t66 = t62 - t4
      t67 = t40 + t12
      t5 = -t3 * t12 * t24 + s16 * t25 + s34 * (((t46 * t30 + ddreal(3) 
     &* t63) * s26 + s25 * (t30 * (t62 - t4) + t36)) * s34 + t7 * (-t3 *
     & t36 * t15 + s26 * (s25 * ((-t22 + t10) * s25 - t60) + s26 * t61))
     &) + (((t3 * s16 * t67 + s34 * t61 - t26 + t18 * (-t5 + t21)) * s56
     & + s25 * (((-t4 + t65) * s25 - ddreal(12) * t12) * s25 + t38) + s3
     &4 * ((s25 * t66 + t9 * s26) * s34 - ddreal(8) * t18 * t15 + t20 * 
     &t26 + t34 * t61) + t39 * t59) * s56 + ddreal(3) * t54 * t59 + t36 
     &* ((s25 * t8 - t65) * s25 + t60) + t4 * t35 * t58 - ((-t8 * t43 * 
     &s16 + t20 * t23 * s26 - s34 * t23 - ddreal(3) * t63) * s34 - s25 *
     & (-t18 * (ddreal(9) * s25 + t9) + (-ddreal(8) * t3 * t64 + ddreal(
     &8) * t23) * s26) - ddreal(3) * t30 * t61 - ddreal(3) * t24) * s34 
     &- t25) * s56 + t13 * t26 + t54 * t58 * s25 + t33 * t59 * s26 * t30
     & + t23 * (t33 * (-s25 + t50) + t41) * s26
      t15 = s16 - t10
      t25 = s16 + t4
      t33 = -s16 + s25 + s34
      t50 = t12 * s25
      t49 = t49 - t12
      t58 = s16 * t1
      t59 = s25 * s34
      t60 = (s56)**2
      t28 = -(t31 * (t40)**2 - ((s26 * (t34 + t21 + t46) * s34 + t3 * (-
     &s26 * t49 + t50) + t27 - ddreal(3) * t47 * t30) * s34 - t58 * t7 *
     & (-t28 - t34 + t46)) * s34) * s56 - t41 * (s34 * t30 + t36 * t3 - 
     &t30 * (s26 + t25)) - ((s26 * t11 + (-s16 - t34 - t4 + s34) * s34 +
     & t49) * s34 + t58 * (s16 - t34 - t21)) * s34 * t60 - t58 * (-t4 + 
     &s16 - s26) * (t7)**2 * s34 + t7 * (t51 + (-t48 + t64) * s26 + t50)
     & * t40 + s34 * (-t46 * s34 + t12 - t18 + t40 - t59) * s56 * t60
      t46 = ddreal(1) / t42
      t45 = ddreal(1) / (t45)**2
      t1 = t59 * (-t23 * t66 + ((s26 + t4) * s34 + (-t62 - t21) * s25 + 
     &(-t25 * t3 + s26) * s26) * s34 + (t33 * s56 - t20 * s25 * (s16 + s
     &34) + t3 * (-s16 * s34 + t33 * s26) + ddreal(3) * t23 + t67) * s56
     & + t50 + t1 * s26 * (-s26 + t15)) * MBL1101(0)
      t4 = (s25 + s26 + s56) * (t37 + (-t15 + s26) * s26 + (s26 - t4) * 
     &s34 + (t10 + t31 - s16 + s34 + s56) * s56 - t18) * MBL1111D2(0)
      t15 = t46 * t45 * ((-t57 + (-t53 - t52) * s56 - t54 - t55 + t56) *
     & MBL1001(0) + s16 * (t57 - (-t53 - t52) * s56 + t54 + t55 - t56) *
     & I300s16(1))
      result = -t45 * (t46 * ((t5 * MBL1011(1) - t28 * MBL1110(1) - t2 *
     & t36 * t44 - ((t42 * s56 + (-ddreal(8) * t12 + t37) * s25 - ((-t22
     & - t34 - t21 + s34) * s34 + t17 * t8 + t14 + t29) * s34 + t36 + t3
     &8 + t39) * s56 + t35 * (t3 * t6 + t34) + t19 * t26 - ddreal(19) * 
     &t12 * t23 + t16 * t13 + ((-t32 + ddreal(3) * t30 + ddreal(3) * t23
     & + ddreal(8) * s25 * (s16 + s26) + ddreal(14) * s16 * s26) * s34 -
     & s25 * ((ddreal(16) * s16 + t10) * s25 + ddreal(17) * t12) - (t3 *
     & (t29 + t14) + t17 * t8) * s26) * s34 + t24) * s56 - (t27 + (s25 *
     & (t19 + t10) + (t22 + t21 + s26) * s26) * s26) * t40 + t43 * t41 -
     & t35 * (s25 * (ddreal(10) * s16 + s25) + (t6 + s26) * s26) + t7 * 
     &((s25 * t11 + t17 * t3 + t14) * s26 + t18 * (t16 + t9)) * s34 + dd
     &real(10) * t12 * t26) / s16 + t1) + t4) / ddreal(2) + t15

           intHLs16s25s26s34s56x1222D6eps0 = result
       end function intHLs16s25s26s34s56x1222D6eps0

       function intHLs16s25s26s34s56x1222D6eps1()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1222D6eps1

           type(dd_complex) :: result

      result = ddreal(0)

           intHLs16s25s26s34s56x1222D6eps1 = result
       end function intHLs16s25s26s34s56x1222D6eps1

       function intHLs16s25s26s34s56x1231D6eps0()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1231D6eps0
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           type(dd_complex) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           type(dd_complex) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           type(dd_complex) ::  t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52
           type(dd_complex) ::  t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63
           type(dd_complex) ::  t64,t65,t66,t67,t68,t69,t7,t8,t9

           type(dd_complex) :: result

      t1 = s25 + s26 + s56
      t2 = s25 - s16
      t3 = ddreal(2)
      t4 = s25 + s26
      t5 = t3 * s26
      t6 = s25 - t5
      t7 = s26 + s16
      t8 = (s25)**2
      t9 = (t8)**2
      t10 = s25 * t8
      t11 = (s26)**2
      t12 = (t11)**2
      t13 = s26 * t11
      t14 = ddreal(3)
      t15 = t3 * s34
      t16 = ddreal(6) * s25
      t17 = s26 * t4
      t18 = t17 * s25
      t19 = t6 * s34
      t20 = s25 + s26 - s34 + s56
      t21 = s16 + s25
      t22 = s16 * s25
      t23 = (-s34 + t21 + s26) * s26 + (t5 + s25 + s16 - s34 + s56) * s5
     &6 + t22
      t24 = ddreal(14) * s16
      t25 = t3 * s16
      t26 = -s16 + s34
      t27 = s26 * t26
      t28 = ddreal(8) * s25
      t29 = ddreal(9) * s26
      t30 = ddreal(4)
      t31 = t30 * s25
      t32 = (s56)**2
      t33 = (t32)**2
      t34 = s56 * t32
      t35 = ddreal(6) * s16
      t36 = t14 * s16
      t37 = t14 * s26
      t38 = t3 * s25
      t39 = ddreal(5) * s16
      t40 = (t38 + t37) * s26
      t41 = ddreal(5) * s26
      t42 = t14 * t11
      t43 = ddreal(5) * s25
      t44 = s16 * s26
      t45 = s16 * t8
      t46 = (s16)**2
      t47 = t46 + t8
      t48 = ddreal(7)
      t49 = t48 * s16
      t50 = -t49 * s25 + t47
      t51 = t14 * t47
      t52 = ddreal(16) * s16
      t53 = -t52 * s25 + t51
      t54 = t30 * s16
      t55 = s25 * t14
      t56 = -t31 + t36
      t57 = ddreal(9) * t46
      t58 = s25 * t47
      t59 = ddreal(12)
      t60 = s26 * t2
      t61 = t14 * t58
      t62 = s16 * t59
      t63 = t50 * t30
      t64 = ddreal(9) * s16
      t65 = (s34)**2
      t66 = s34 * t65
      t67 = ddreal(10) * s16
      t68 = t3 * s34 * t2 - (s25 * t48 + t41) * s16 + t47
      t69 = t4 * s34
      t51 = t69 * (((-t56 * s26 - s25 * (-t54 + t55)) * s26 + t45) * s34
     & + t4 * ((-t28 * s16 + t51) * s26 - t3 * (t11 * t2 + t45)))
      t6 = s25 * (-t50 * t8 + t6 * t66) * s26
      t6 = s16 * s25 * t9 + s16 * s26 * t12 + s16 * s56 * t33 - ((-t48 *
     & t9 - ddreal(5) * t12) * s16 + ((-s25 * (-s26 * t30 + s25) * s34 +
     & t11 * t14 * t56 + ((-t39 + t55) * s25 - t5 * t21) * s25) * s34 - 
     &t4 * ((-ddreal(8) * t60 + (-ddreal(20) * s16 + t43) * s25 + t57) *
     & s26 + t61 - t62 * t8)) * s34 + t10 * t47 + t63 * t13 + t40 * s25 
     &* t53) * s56 - t12 * t50 - t32 * (t14 * (s26 * t53 + t58) * s25 + 
     &s34 * (-(s25 * t21 + (s25 * t59 - t64) * s26) * s34 + s25 * ((t52 
     &- t31) * s25 - ddreal(6) * t46) + (t60 * t59 + (t62 + t55) * s25 -
     & t57) * s26) - t52 * t10 - t67 * t13 + ddreal(6) * t50 * t11 + t38
     & * t66) - t33 * t68 - t34 * ((-t67 * s26 + t63) * s26 + t61 + s34 
     &* (s34 * t56 + (s25 + t54) * s25 - t14 * t46 + ddreal(8) * t60) - 
     &t52 * t8) - s25 * t53 * t13 - t8 * t53 * t11 + t51 + t6
      t12 = ddreal(6) * s26
      t21 = t59 * s26
      t47 = (t4)**2
      t50 = s16 * t4
      t51 = s16 * t2
      t52 = ddreal(1) / s16
      t23 = ddreal(1) / (t23)**2
      t20 = ddreal(1) / (t20)**2
      t6 = t6 * MBL1011(1) + ((t25 * (t47)**2 + t42 * t66) * s34 + ((ddr
     &eal(8) * s16 * t4 * t47 + t12 * t66) * s34 - t65 * (((t21 + t49) *
     & s25 + ddreal(18) * s26 * t7) * s34 - t4 * ((ddreal(21) * s16 + t1
     &6) * s26 + t59 * (t11 + t46) - t22))) * s56 + t32 * ((((ddreal(10)
     & * s25 + ddreal(21) * s26 + t35) * s16 + t14 * t8 + ddreal(18) * t
     &17) * s34 + t62 * t47) * s34 + t14 * t66 * (s34 - t12 - t36 - t38)
     &) + s34 * ((-ddreal(6) * s34 + t21 + t49 + t16) * s34 + ddreal(8) 
     &* t50) * t34 + t4 * ((-t12 - t64) * s26 + t38 * s16) * t66 - t47 *
     & ((-s26 * t48 + t31 - t35) * s16 - t42) * t65 + s34 * (s34 * t14 +
     & t25) * t33) * MBL1110(1)
      t4 = t20 * ((t10 * t3 - ((s16 - t31 - s26) * s26 - s25 * (-t25 + t
     &43)) * s26 - (-(t38 - t37) * s34 - t4 * (-t31 + t36 + t5)) * s34 -
     & ((s16 - t15 - t31 - t37 - s56) * s56 - (s34 * t30 + t28) * s26 - 
     &t14 * (-s34 * t26 + t11) + t3 * ((s16 + s34) * s25 + t44) - ddreal
     &(5) * t8) * s56 - t45) * MBL1001(0) + t14 * t18 + ((-t2 * t3 + s26
     &) * t4 + t19) * s34 + ((t14 * t4 + s34 + s56) * s56 + (t3 * t7 - s
     &25 - t15) * s34 + t14 * (t11 + t8) + t16 * s26) * s56 + t10 + t13)
     & * t1
      t1 = ddreal(1) / t1
      t1 = t52 * (((t14 * ((-(-t5 * s16 * s34 - ((s26 + s34 + s25) * s25
     & - t50) * s16 + s26 * t65) * s56 + t3 * t44 * (s25 * t2 + t69) + s
     &16 * t10 - t11 * (-t51 + t65) - t46 * t8) * s56 + t51 * t18) + s16
     & * t9 - t10 * t46 - t13 * t65 - t34 * (-t25 * s34 - t51 + t65) + t
     &51 * t13 - t19 * s16 * t47) * MBL1111D2(0) + s25 * s34 * (t3 * ((-
     &t27 + t22) * s56 + t22 * s26) - t11 * t26 - t26 * t32 + t45) * MBL
     &1101(0)) * t1 * t23 + t65 * MBL1110(0) * t20)
      result = t23 * (-((t41 * s25 + t3 * t8) * s16 + (t39 * s25 + (t14 
     &* (s26 + s16 - s34) + s25 + s56) * s56 - ddreal(6) * t27 + t40) * 
     &s56 + t11 * (s25 + t36 + s26) - t42 * s34) * I300s16(1) + (t6 * t2
     &0 - t14 * (t34 + t13) - (s25 * (s25 + t24) - s34 * (s25 + ddreal(8
     &) * s26)) * s26 - ((-ddreal(8) * t26 + t31 + t29) * s56 + (t29 + t
     &28) * s26 - ddreal(16) * t27 + s25 * (s25 - s34) + t24 * s25) * s5
     &6 - t35 * t8 - t30 * (s25 + t25) * t11) * t52 + t4) / ddreal(4) + 
     &t1 / ddreal(2)

           intHLs16s25s26s34s56x1231D6eps0 = result
       end function intHLs16s25s26s34s56x1231D6eps0

       function intHLs16s25s26s34s56x1231D6eps1()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1231D6eps1
           type(dd_complex) ::  t1,t2

           type(dd_complex) :: result

      t1 = s25 + s26 - s34 + s56
      t2 = -ddreal(1) / ddreal(4)
      t1 = ddreal(1) / (t1)**2
      result = ddreal(1) / s16 * t1 * (t2 * (s25 + s26 - ddreal(3) * s34
     & + s56) + (s34)**2 * MBL1110(1) / ddreal(2))

           intHLs16s25s26s34s56x1231D6eps1 = result
       end function intHLs16s25s26s34s56x1231D6eps1

       function intHLs16s25s26s34s56x1311D4eps0()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1311D4eps0
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           type(dd_complex) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           type(dd_complex) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           type(dd_complex) ::  t5,t6,t7,t8,t9

           type(dd_complex) :: result

      t1 = (s16)**2
      t2 = s16 * t1
      t3 = (s25)**2
      t4 = s16 + s25
      t5 = ddreal(2)
      t6 = s16 * s25
      t7 = s34 * t4
      t8 = (s34)**2
      t9 = t3 + t1 + t8
      t10 = s26 * s34
      t11 = t1 * s25
      t12 = s16 * t3
      t13 = (t3 + t1) * s26
      t14 = (-s26 * t4 * t5 + t10 - t6) * s34 + (-t5 * t7 + t9) * s56 + 
     &t11 + t12 + t13
      t7 = -t5 * (t7 + t6) + t9
      t9 = t5 * s26
      t15 = (-s34 + t4 + s26) * s26 + (t9 + s16 + s25 - s34 + s56) * s56
     & + t6
      t16 = t5 * s25
      t17 = ddreal(3) * s16
      t18 = t5 * s16
      t19 = s25 - t18
      t20 = ddreal(4)
      t21 = (s56)**2
      t22 = (s26)**2
      t23 = s25 + s26 - s34 + s56
      t24 = s26 + s56
      t25 = t4 * t3
      t26 = ddreal(3) * t1
      t27 = t26 * s25
      t28 = s16 - s25
      t29 = s25 * t28
      t30 = s25 + s26
      t31 = t4 * t8
      t32 = s26 * t28
      t33 = s16 + t16
      t34 = s25 + t18
      t35 = ddreal(5) * s25
      t36 = t34 * s26
      t37 = s16 * (t35 - t18)
      t38 = (ddreal(7) * s16 + s25) * s25
      t39 = -ddreal(9)
      t40 = ddreal(6) * s16
      t23 = ddreal(1) / t23
      t15 = ddreal(1) / t15
      t7 = ddreal(1) / t7
      t41 = ddreal(1) / s16
      t13 = t24 * (t12 * (-t17 + t16) + s25 * t2 + (t31 - t27 + t25 + t2
     &) * s26 - (t5 * (t13 + t12) + s25 * (t22 + t1)) * s34 + (-t5 * ((s
     &25 * t30 + t1) * s34 + t32 * s25) + t2 + t25 - t27 + t31) * s56 - 
     &t22 * t29 - s25 * (s16 - s25 + s34) * t21) * MBL1011(1)
      t25 = s25 * ((t32 + t1) * s25 - s26 * t8 + (s26 * t33 + t6) * s34 
     &+ ((t33 - s34) * s34 + t29) * s56 - t12) * MBL1101(0)
      t10 = (-t5 * s25 * (t22 + t21) + ((-t17 - t16) * s25 + t1) * s26 +
     & (s26 * t19 + t10 - t6) * s34 + ((-s26 * t20 - t16 - t17) * s25 + 
     &(t19 + s34) * s34 + t1) * s56 + t11 - t17 * t3) * MBL1001(0) * t23
      t19 = (s25 + s26 + s56) * MBL1110(0) * t41 * t23
      t1 = t15 * (t7 * (t14 * I300s16(1) + (t23 * ((-t1 * s25 * t3 + (((
     &-s25 * s34 - t18 * s34 + t1 - t6 + t8) * s56 + (-t3 * t5 + ddreal(
     &3) * t32) * s16 + (-ddreal(3) * s25 * s26 - t40 * s26 + (ddreal(3)
     & * s16 + ddreal(3) * s26 + t16 - s34) * s34 - ddreal(3) * t1 - t38
     &) * s34 + t1 * t4) * s56 + t5 * s16 * (t32 * t33 + t11) - ((-(s16 
     &- t9) * s34 - (s25 * t20 + ddreal(3) * s26 + t40) * s26 - t37) * s
     &34 - s16 * ((s16 * t39 - t35) * s25 + t1) - (-t5 * (t26 + t38) - d
     &dreal(3) * t36) * s26) * s34 - t12 * t4 + t17 * t28 * t22) * s56 +
     & t2 * t3 - t8 * (((-t17 - t16 - s26) * s26 - t37) * s26 + t11) - t
     &30 * (ddreal(6) * t11 - (-s16 * (t35 + t17) - t36) * s26 - t2) * s
     &34 + s26 * (s16 - s26) * s34 * t8 + t32 * (s25 * t34 + (t33 + s26)
     & * s26) * s16) * MBL1110(1) + t13) + t25) * t41 + t10) - (t24)**2 
     &* MBL1111D2(0) * t41) + t19
      result = t1 / ddreal(2) + t14 * t41 * t7 * t15

           intHLs16s25s26s34s56x1311D4eps0 = result
       end function intHLs16s25s26s34s56x1311D4eps0

       function intHLs16s25s26s34s56x1311D4eps1()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1311D4eps1
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s25 + s26 - s34 + s56
      t1 = ddreal(1) / t1
      result = ddreal(1) / s16 * t1 * ((s25 + s26 + s56) * MBL1110(1) + 
     &ddreal(1)) / ddreal(2)

           intHLs16s25s26s34s56x1311D4eps1 = result
       end function intHLs16s25s26s34s56x1311D4eps1

       function intHLs16s25s26s34s56x1312D6eps0()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1312D6eps0
           type(dd_complex) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           type(dd_complex) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           type(dd_complex) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           type(dd_complex) ::  t131,t14,t15,t16,t17,t18,t19,t2,t20,t21,t22,t23
           type(dd_complex) ::  t24,t25,t26,t27,t28,t29,t3,t30,t31,t32,t33,t34
           type(dd_complex) ::  t35,t36,t37,t38,t39,t4,t40,t41,t42,t43,t44,t45
           type(dd_complex) ::  t46,t47,t48,t49,t5,t50,t51,t52,t53,t54,t55,t56
           type(dd_complex) ::  t57,t58,t59,t6,t60,t61,t62,t63,t64,t65,t66,t67
           type(dd_complex) ::  t68,t69,t7,t70,t71,t72,t73,t74,t75,t76,t77,t78
           type(dd_complex) ::  t79,t8,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89
           type(dd_complex) ::  t9,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99

           type(dd_complex) :: result

      t1 = s16 - s25
      t2 = ddreal(2)
      t3 = (s16)**2
      t4 = (t3)**2
      t5 = s16 * t4
      t6 = s16 * t3
      t7 = t2 * t3
      t8 = (s25)**2
      t9 = (t8)**2
      t10 = t8 * t9
      t11 = s25 * t9
      t12 = s25 * t8
      t13 = t2 * t8
      t14 = ddreal(3)
      t15 = t14 * s25
      t16 = s25 + s26
      t17 = -ddreal(37) * s16
      t18 = ddreal(11) * s25
      t19 = ddreal(8) * s16
      t20 = ddreal(7) * t3
      t21 = t2 * t6
      t22 = ((-t19 + s25) * s25 - t20) * s25 - t21
      t23 = t3 * t12
      t24 = ddreal(7) * s16
      t25 = t6 * s25
      t26 = ddreal(15)
      t27 = ddreal(17) * s16
      t28 = ddreal(4)
      t29 = ddreal(5) * s16
      t30 = t28 * s25
      t31 = ddreal(11) * t3
      t32 = (t29 - t30) * s25
      t33 = -ddreal(22) * t6
      t34 = (t32 + t31) * s25 + t33
      t35 = t14 * s16
      t36 = t35 - s25
      t37 = (s26)**2
      t38 = s26 * t37
      t39 = s16 * s25
      t40 = t28 * s16
      t41 = ddreal(18)
      t42 = ddreal(6)
      t43 = s25 * t42
      t44 = t41 * t3
      t45 = (t29 + t43) * s25 + t44
      t46 = t2 * s16
      t47 = t46 + t15
      t48 = t3 * t8
      t49 = t14 * s26
      t50 = t39 * t2
      t51 = -ddreal(44) * t3
      t52 = ddreal(21) * t6
      t53 = ddreal(13) * t4
      t54 = ddreal(20)
      t55 = ddreal(9) * s25
      t56 = s16 * t26
      t57 = ddreal(12) * s16
      t58 = ddreal(14) * s16
      t59 = ddreal(8) * s25
      t60 = (t1)**2
      t61 = (t60)**2
      t62 = t4 * s25
      t63 = t60 * (s16 * t1 + t13)
      t64 = t6 * t8
      t65 = (-t15 + s16) * t60
      t66 = (s34)**2
      t67 = (t66)**2
      t68 = s34 * t67
      t69 = s34 * t66
      t70 = t3 * t9
      t71 = t2 * s26
      t72 = s16 * t12
      t73 = ddreal(9) * s26
      t74 = ddreal(11) * s16
      t75 = s16 * s26
      t76 = ddreal(9) * t6
      t77 = s16 * t8
      t78 = ddreal(7) * s25
      t79 = t3 * s25
      t33 = (t14 * s25 * (t36 * t66 - t77) + (t66 * (t47 - s34) + t22) *
     & s34 + t4 + t79 * (-t29 + t78)) * s56 - t26 * t48 * t1 - t14 * s16
     & * (-t65 * s26 + t4) + ((((-t49 - t24 - t30 + s34) * s34 + (t73 + 
     &t29) * s25 + t42 * (t75 + t8) + t44) * s34 + (s16 * (ddreal(27) * 
     &s26 + t74) + (-t73 + t29 - t30) * s25) * s25 + t33) * s34 + ((t51 
     &+ t8) * s25 - t52) * s25 + t14 * (t22 * s26 + t72) + t53) * s34 + 
     &t39 * (-t12 * t42 + t76)
      t44 = t16 * (t25 * (-ddreal(22) * s25 + t24) - (-s16 * (((t17 + t1
     &8) * s25 - ddreal(19) * t3) * s25 + ddreal(13) * t6) - t22 * s26) 
     &* s26 - t23) * s34
      t18 = (t33 * s56 - t14 * s16 * (-t65 * t37 + t11) + t26 * t64 * t1
     & - t42 * s16 * (t63 * s26 + t62) + ((t2 * t34 * s26 + t55 * t36 * 
     &t37 + (-((t58 + t49 + t59) * s26 + t50) * s34 + (t14 * t47 * s26 +
     & t2 * t45) * s26 + t39 * (t55 + t57)) * s34 - t56 * t12 - ddreal(2
     &4) * t25 + ddreal(17) * t48) * s34 + (t2 * ((((t35 + s25) * s25 + 
     &t51) * s25 - t52) * s25 + t53) + t49 * t22) * s26 + t39 * (((-ddre
     &al(38) * s16 + t18) * s25 - ddreal(41) * t3) * s25 + t54 * t6)) * 
     &s34 + ddreal(9) * t70 + t71 * t68) * s56 + t11 * t3 + t37 * t68 + 
     &((t38 * t36 * t14 - ddreal(5) * t25) * s25 + t34 * t37 - t23 - t39
     & * ((s25 * t26 - t27) * s25 + ddreal(24) * t3) * s26) * t66 + (t49
     & * t39 * (t40 + t15) + t37 * (s26 * t47 + t45) + t48) * t69 + t44 
     &+ t65 * s16 * t38 - s26 * ((t24 + t30 + s26) * s26 + t50) * t67
      t22 = t4 * t8
      t33 = t60 * s25
      t34 = ddreal(5) * s25
      t44 = t6 * t12
      t45 = s16 + s25
      t51 = -t2 * (s34 * t45 + t39) + t66 + t8 + t3
      t52 = (-s34 + t45 + s26) * s26 + (t71 + s16 + s25 - s34 + s56) * s
     &56 + t39
      t53 = t8 + t3
      t65 = t2 * t53
      t68 = t2 * s25
      t80 = -t68 + s16
      t81 = ddreal(8) * t3
      t82 = t28 * t6
      t50 = t14 * t53 + t50
      t83 = -ddreal(10) * s25
      t84 = -ddreal(10) * s16
      t85 = t42 * t3
      t86 = (t84 - t34) * s25 + t85
      t87 = t68 + s16
      t88 = t14 * t87
      t89 = t2 * s25 * t45
      t90 = (t40 - t15 - s26) * s26 + t89
      t91 = ddreal(8) * t6
      t92 = ((-ddreal(40) * s16 + t43) * s25 - ddreal(22) * t3) * s25 + 
     &t91
      t93 = t41 * t77 * t45
      t94 = t42 * s25 * (t12 + t6)
      t95 = s25 * t53
      t96 = (t2 * t86 - t73 * t87) * s26 + t42 * t95 - ddreal(11) * t77
      t89 = (t19 - t43 - t49) * s26 + t89
      t97 = s16 * s34
      t98 = t60 * ((-t46 + t30) * s25 + t3)
      t99 = t48 * t53
      t100 = t85 * t12
      t101 = t14 * t80 * t60 * t37
      t102 = t14 * t50 * s26
      t103 = s16 * t42
      t104 = s25 * (-t12 + t6)
      t105 = ddreal(9) * s16
      t54 = -t104 * t28 + (((t14 * t16 + s34 - t40) * s34 + s16 * (-t73 
     &+ t103) + (-s26 * t41 - t34 + t84) * s25) * s34 + ((s16 * t54 - t1
     &5) * s25 + t31) * s25 - t82 + t102) * s34 + t4 - t49 * t80 * t60 +
     & t77 * (t105 + t83)
      t106 = t28 * t3
      t107 = s25 * ((t29 - t68) * s25 - t106) + (-t66 - t50) * s34 + t6 
     &+ t88 * t66
      t108 = t98 + t67
      t109 = (s56)**2
      t110 = s56 * t109
      t83 = ((-t88 * t38 + (s25 * (t42 * t53 - ddreal(11) * t39) + t86 *
     & s26) * s26 + t48) * s34 - t16 * (t79 * (t83 + t46) + s26 * (-t50 
     &* s26 + ((-s16 * t41 + t43) * s25 - t81) * s25 + t82))) * s34
      t86 = t33 * (t65 - t39)
      t88 = t108 * t37
      t111 = s26 + s56
      t36 = t36 * t87
      t112 = ddreal(16) * t79
      t113 = t42 * (t1 * t8 + t6)
      t114 = t113 - t112
      t115 = t14 * t3
      t116 = (t103 + s25) * s25 + t115
      t117 = (t29 + t68) * s25
      t118 = t117 - t106
      t119 = t75 - t8
      t120 = ddreal(12) * s25
      t121 = ddreal(26) * t3
      t122 = t6 * t42
      t123 = t36 * t49
      t124 = t49 * t116
      t125 = t5 * s25
      t126 = ddreal(16) * t3
      t127 = ddreal(13) * t3
      t128 = t4 * t42
      t13 = (t28 * t39 * t53 + ((t97 - t116) * s34 + t36 * t1) * s34 - t
     &4 - t9 - t48 * t42) * s56 + s25 * ((((t105 - t68) * s25 - t126) * 
     &s25 + ddreal(14) * t6) * s25 - t128) - (((-t97 - (t34 + t49 - t40)
     & * s16 - t13) * s34 + t112 - t113 + t124) * s34 + t1 * (((t43 - t5
     &6) * s25 - t127) * s25 - t123 + t82)) * s34 + t5 - t49 * t61
      t40 = t5 * t8
      t112 = t1 * t16
      t36 = t112 * (t2 * (t45 * (-ddreal(7) * t39 + t65) * s26 + t25) - 
     &t36 * t37 - t77 * (t29 + t15)) * s34
      t10 = s16 * t10 + (t13 * s56 + t2 * (t80 * t61 * s26 + t125) - (((
     &-t2 * t118 * s26 - (t71 * s16 - t8) * s34 - t12 * t28 + t79 * t2 -
     & t35 * t37 + ddreal(5) * t77) * s34 - (t28 * t114 / ddreal(2) - t1
     &24) * s26 + t42 * s25 * ((t8 + t3) * s25 - t6) - ddreal(16) * t72)
     & * s34 + t1 * (s25 * (((-ddreal(13) * s16 + t30) * s25 - t26 * t3)
     & * s25 + t122) + (-t123 + ((-ddreal(30) * s16 + t120) * s25 - t121
     &) * s25 + t91) * s26)) * s34 - t10 - t14 * t61 * t37 + t77 * (((-t
     &58 + t43) * s25 + t126) * s25 - t76)) * s56 + t66 * ((s26 * t119 *
     & s34 + (s25 * (-t32 - t7) + (t75 + t118) * s26) * s26 - t72) * s34
     & + s25 * (t77 * (t103 + t15) + t71 * (-t14 * (t95 - t6) + ddreal(8
     &) * t77)) + t37 * (-s26 * t116 + t114) + t64) + t40 - t36 + (s25 *
     & (t46 - s25) + (t80 - s26) * s26) * s26 * t61
      t13 = t39 * t28 + t53
      t32 = (-t78 + s16) * s25 + t115
      t36 = t55 * t37
      t45 = t28 * s26
      t61 = t14 * t8
      t65 = s26 * s34
      t56 = (((t68 + t56) * s25 - t106) * s25 - t91) * s25
      t66 = (t34 - t46) * s25 + t115
      t72 = t66 * s26
      t76 = ((t43 + t24) * s25 - ddreal(5) * t3) * s25 + t14 * t6
      t91 = t3 + ddreal(5) * t8
      t95 = t30 + t35
      t113 = (t59 + s16) * s25
      t114 = t113 + t106
      t116 = t78 + t35
      t118 = ((ddreal(13) * s25 + t105) * s25 - t31) * s25
      t41 = t41 * s25
      t123 = (t118 + t6) * t14
      t124 = t45 * t66
      t126 = t95 * t37
      t117 = t1 * (-t117 + t3)
      t129 = t1 * t13
      t130 = -ddreal(29)
      t131 = t117 * t37
      t26 = ((s25 * ((t68 + t35) * s25 - t85) + ((t95 - s34) * s34 - t66
     &) * s34 + t6) * s56 + t28 * (t117 * s26 - t104) - ddreal(12) * t77
     & * t1 + ((ddreal(16) * s25 * s26 + (-t78 - t45 - t35 + s34) * s34 
     &+ t26 * t8 + t115 + ddreal(12) * t75) * s34 - t118 - t124 - t6) * 
     &s34) * s56 + s25 * (-ddreal(12) * t129 * s26 + t56) + t42 * (-t62 
     &+ t131) + ((((t49 + t68 + s16) * s34 - (s26 * t42 + ddreal(21) * s
     &25 + t105) * s26 - t106 - t113) * s34 + s25 * ((t120 + t58) * s25 
     &- ddreal(10) * t3) + t42 * (t126 + t6) + t73 * t91) * s34 + (((s16
     & * t130 - t59) * s25 + t121) * s25 + t26 * t6) * s25 + (-t72 * t42
     & - t123) * s26 - t28 * t4) * s34 + t5
      t11 = t26 * s56 + t2 * (t125 + t1 * ((((-t27 - t68) * s25 - t127) 
     &* s25 - ddreal(5) * t6) * s25 + t4) * s26) + t28 * (t117 * t38 + t
     &70) + (((-t2 * (s26 * t114 + t79) - t37 * (t14 * t116 + t45) - t77
     & * t42 + t65 * (t46 + t49 + t30)) * s34 + t28 * s26 * (t126 + t76)
     & + ddreal(9) * t91 * t37 + t39 * ((t41 + t84) * s25 + t85)) * s34 
     &+ ((-t124 - t123) * s26 - t1 * (((t17 - t59) * s25 - t31) * s25 + 
     &t82) * t2) * s26 + t39 * (((t19 - t41) * s25 + ddreal(28) * t3) * 
     &s25 - t122)) * s34 - ddreal(12) * s25 * (t129 * t37 + t62) + t103 
     &* t11
      t17 = s25 + s26 + s56
      t19 = ddreal(1) / (t52)**2
      t26 = ddreal(1) / (t51)**2
      t5 = s34 * (s25 * (-t129 * t28 * t38 + t64 * (-t103 + t15)) + t2 *
     & t39 * (s16 * t9 + ((t47 * t8 - t122) * s25 + t4) * s26) + (((-s26
     & * (((t116 + s26) * s26 + t114) * s26 + t39 * (t43 + t46)) + t37 *
     & (t68 + s16 + s26) * s34) * s34 + t2 * (t37 * t76 + t23 + t39 * ((
     &-t29 + t55) * s25 + t115) * s26) + (t37)**2 * t95 + t64 + t14 * t9
     &1 * t38) * s34 - t16 * (t2 * s16 * ((((t55 - t103) * s25 - t81) * 
     &s25 + t21) * s26 + t25) + t37 * (((t74 + t59) * s25 - ddreal(14) *
     & t3) * s25 + t6 + t72) + t48 * (-t57 + t30))) * s34 + t11 * s56 + 
     &t37 * (t131 + (t56 - t128) * s25 + t5) + t40) * MBL1110(1) - t14 *
     & s16 * (t63 * t37 + t22 + t33 * (-s25 * t1 + t7) * s26) + t18 + t4
     &4 * (-t34 + t24) + t111 * (-t28 * t23 * t53 + t122 * t9 + t10) * M
     &BL1011(1)
      t1 = t19 * ((-t61 * (s16 - s25 + s34) * s34 * t110 + s25 * s34 * (
     &s34 * (t65 * (-(s25 - s26) * s34 - t119 * t14 - t2 * s25 * (s16 - 
     &s26)) + t14 * s25 * (t79 - t38) - s26 * (s25 * ((t15 + s16) * s25 
     &- t20) - t32 * s26)) + (((-t73 * t1 + t28 * t8 - t115) * s25 - ((-
     &t68 + t35 - s34) * s34 - (-t73 - t78 + s16) * s25 - t115) * s34 - 
     &t6) * s56 - t36 * t1 - t71 * t1 * (t87)**2 + t39 * (-t106 + t61) +
     & (((t71 - s25) * s34 + (-t46 + t45 + t15) * s25 - t75 * t42) * s34
     & - t12 * t14 + t71 * t32 - t36 - t77 + ddreal(7) * t79) * s34 + t9
     &) * s56 - t112 * (s26 * t13 + t14 * s25 * (t37 + t3)))) * MBL1101(
     &0) * t26 + t111 * (t17)**2 * MBL1111D2(0))
      result = t26 * t19 * ((t2 * t4 * t12 - s16 * (t83 + t99) + ((t101 
     &- t100) * s16 - t2 * s16 * ((t9 + t4) * s25 + (t98 + t67) * s26) +
     & t97 * ((s34 * t89 - t96) * s34 + (-t49 * t50 + t92) * s26 - t93 +
     & t94) + ddreal(5) * t99) * s56 + s16 * t80 * t60 * t38 + s16 * t10
     &7 * t110 - s16 * t54 * t109 - t88 * s16 - t75 * (-t69 * t90 + t86)
     &) * I300s16(1) + (-t44 * t2 + (-t69 * t90 + t86) * s26 + ((-t107 *
     & s56 + t54) * s56 + t2 * ((t9 + t4) * s25 + s26 * t108) + t100 - t
     &101 - ddreal(5) * t77 * t53 + ((-s34 * t89 + t96) * s34 + t93 - t9
     &4 + s26 * (t102 - t92)) * s34) * s56 + t22 + t83 + t88 + t70 - t80
     & * t60 * t38) * MBL1001(0) + t5 / s16) / ddreal(4) + t1 / ddreal(2
     &)

           intHLs16s25s26s34s56x1312D6eps0 = result
       end function intHLs16s25s26s34s56x1312D6eps0

       function intHLs16s25s26s34s56x1312D6eps1()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1312D6eps1

           type(dd_complex) :: result

      result = ddreal(0)

           intHLs16s25s26s34s56x1312D6eps1 = result
       end function intHLs16s25s26s34s56x1312D6eps1

       function intHLs16s25s26s34s56x1321D6eps0()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1321D6eps0
           type(dd_complex) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           type(dd_complex) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           type(dd_complex) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           type(dd_complex) ::  t131,t132,t133,t134,t135,t136,t137,t14,t15,t16,t17,t18
           type(dd_complex) ::  t19,t2,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29
           type(dd_complex) ::  t3,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t4
           type(dd_complex) ::  t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t5,t50
           type(dd_complex) ::  t51,t52,t53,t54,t55,t56,t57,t58,t59,t6,t60,t61
           type(dd_complex) ::  t62,t63,t64,t65,t66,t67,t68,t69,t7,t70,t71,t72
           type(dd_complex) ::  t73,t74,t75,t76,t77,t78,t79,t8,t80,t81,t82,t83
           type(dd_complex) ::  t84,t85,t86,t87,t88,t89,t9,t90,t91,t92,t93,t94
           type(dd_complex) ::  t95,t96,t97,t98,t99

           type(dd_complex) :: result

      t1 = s16 - s25
      t2 = s25 * t1
      t3 = ddreal(3)
      t4 = ddreal(2)
      t5 = (s16)**2
      t6 = (t5)**2
      t7 = s16 * t5
      t8 = t3 * t5
      t9 = t4 * t2
      t10 = t8 - t9
      t11 = t4 * s16
      t12 = s25 + t11
      t13 = s16 + s25
      t14 = s25 * t13
      t15 = ddreal(9) * t5
      t16 = -t15 + t14
      t17 = ddreal(6)
      t18 = (s26)**2
      t19 = (t18)**2
      t20 = s26 * t18
      t21 = (s25)**2
      t22 = (t21)**2
      t23 = s25 * t22
      t24 = s25 * t21
      t25 = ddreal(4) * t13
      t26 = t17 * s16
      t27 = t26 * s25
      t28 = ddreal(4) * s25
      t29 = t4 * s26
      t30 = ddreal(9) * s16
      t31 = t13 * t18
      t32 = t5 * s25
      t33 = ddreal(12)
      t34 = t26 * t21
      t35 = s16 * s25
      t36 = (s34)**2
      t37 = (t36)**2
      t38 = s34 * t36
      t39 = (t1)**2
      t40 = t32 * t1
      t41 = s26 * t1
      t42 = t33 * s26
      t43 = t17 * s26
      t44 = ddreal(4) * s16
      t45 = t4 * s25
      t46 = ddreal(5) * t5
      t47 = t43 * t39
      t48 = s34 * t13
      t49 = -t4 * (t35 + t48) + t36 + t5 + t21
      t50 = t18 * t38
      t51 = t5 * t21
      t52 = (s56)**2
      t53 = t39 * t20
      t54 = s26 * s34
      t55 = (s34 - t13 - s26) * s26 + (-t29 - s16 - s25 + s34 - s56) * s
     &56 - t35
      t56 = t3 * s16
      t57 = -s25 + t56
      t58 = -t5 + t21
      t59 = t11 * s25
      t60 = t59 - t58
      t61 = t3 * s25
      t62 = (t26 - t61) * s25 + t5
      t63 = -s25 + t11
      t64 = s25 + s26
      t65 = t5 + t21
      t66 = ddreal(8) * s16
      t67 = ddreal(5) * s25
      t68 = t67 + t44
      t69 = t68 * s26
      t70 = t32 * t17
      t71 = ddreal(32) * s16
      t72 = ddreal(19) * s25
      t73 = t4 * t7
      t74 = -ddreal(33)
      t75 = ddreal(4) * t5
      t76 = ddreal(8) * s26
      t77 = t3 * t18
      t78 = t1 * t62
      t79 = s16 * t21
      t80 = t18 * t63
      t81 = t17 * s25
      t82 = ddreal(13) * s25
      t83 = ddreal(11) * s25
      t84 = t5 * t33
      t85 = t36 + t5
      t86 = ddreal(5) * s34
      t87 = (t64)**2
      t62 = (((s25 - t86 - t56) * s25 + t4 * t85 - t44 * s34) * s56 - s2
     &5 * ((-t61 + t30) * s25 - t46) + ((t56 + t83 + t76 - s34) * s34 - 
     &ddreal(16) * s16 * t64 - (ddreal(20) * s26 + t82) * s25 - t8) * s3
     &4 + t7 + ddreal(4) * t41 * t63) * s56 + t3 * (s25 * (s25 * t65 + t
     &7) + t36 * ((s16 - s25 - s26) * s34 - (-t56 - t83) * s26 + ddreal(
     &4) * t14 + ddreal(4) * t18 - t5) + t41 * t62) - s34 * (t38 + ((ddr
     &eal(20) * s16 + t83) * s25 + t84) * s25 + (t43 * t68 + t13 * (t56 
     &+ t82) * t3) * s26 - t7) + t80 * t17 * t1 - t30 * t24
      t68 = t18 * t37
      t88 = t5 * t24
      t89 = t61 * t1
      t90 = s25 + s26 - s34 + s56
      t91 = -t5 + t2
      t92 = -t8 + t14
      t93 = s25 + t56
      t94 = t43 * t13
      t95 = t3 * s26
      t96 = t1 * t91
      t97 = s16 - s26
      t98 = t7 * t21
      t92 = t29 * t40 + (t18 * t92 - t4 * s26 * (t12 * t35 - t31) - t51)
     & * s34 + ((-s56 * t49 - t4 * t35 * t1 + ((t3 * t97 + s25 - s34) * 
     &s34 + t14 - t8 + t94) * s34 - t24 + t7 - t95 * t39) * s56 - t4 * (
     &(t96 + t38) * s26 - t40) + s34 * (((t26 - t95) * s26 + t4 * s25 * 
     &(s16 + s26)) * s34 + (t4 * t92 + t94) * s26 - t35 * (t44 + t45)) -
     & t77 * t39) * s56 - t50 - t53 - t88 + t98 + s26 * (t59 - s26 * (s2
     &6 - t93)) * t36 - t96 * t18
      t96 = ddreal(5) * s16
      t99 = t3 * t7
      t100 = t56 - t45
      t101 = s25 * t100 + t5
      t102 = (-t81 + t30) * s25 + t5
      t103 = ddreal(7)
      t104 = t103 * s25
      t105 = t104 + t26
      t106 = t105 * s26
      t107 = ddreal(11) * s16
      t108 = (t107 + t81) * s25
      t109 = ddreal(8) * s25
      t110 = ddreal(8) * t5
      t111 = ddreal(16) * s25
      t112 = t1 * t102
      t113 = ddreal(4) * t100
      t114 = ddreal(5) * t14
      t115 = ddreal(22) * s16
      t116 = ddreal(14) * t5
      t117 = -ddreal(48)
      t118 = ddreal(15) * s16
      t119 = ddreal(18) * s26
      t120 = ddreal(16) * s16
      t121 = t100 * t18
      t122 = t7 * s25
      t123 = t26 * s34
      t85 = t3 * t85
      t102 = (-((s34 * t103 - t45 + t96) * s25 + t123 - t85) * s56 - s25
     & * ((-t81 + t118) * s25 - t110) - t38 * t4 + t36 * (t96 + t111 + t
     &42) + t7 - ddreal(4) * ((t104 + t26) * s26 + t114 + t5) * s34 + dd
     &real(4) * t41 * t100) * s56 + t17 * (t1 * t121 + t21 * t65) + t3 *
     & (t102 * t41 + t122) - (((t43 + t67 - t56 + s34) * s34 + (-t120 - 
     &t72) * s25 + (s25 * t117 - t118 - t119) * s26 + t8) * s34 + ((t115
     & + t72) * s25 + t116) * s25 + (t33 * (t5 + t114) + t43 * t105) * s
     &26 - t7) * s34 - t118 * t24
      t104 = t7 * t24
      t9 = -t89 * t101 * t18 - (-t38 * ((t67 - t56) * s26 - t4 * (-t18 +
     & t2)) - t21 * (t21 * (-t45 + t96) - t99)) * s26 - (t102 * s56 + t3
     & * (t112 * t18 + t98) - t4 * ((s26 * t38 + t64 * (((t44 + t61) * s
     &25 + t110) * s25 + t4 * s26 * (t106 + (t30 + t109) * s25 + t8) - t
     &7)) * s34 - t23) - t36 * ((-s26 * t17 * t97 + ddreal(10) * s25 * s
     &26 - t9) * s34 + s25 * (-t108 + t75) + ((-t3 * (t96 + t111) - t42)
     & * s26 + t13 * (t56 - t72) * t4) * s26) + t81 * t41 * t101 + t113 
     &* t1 * t20 - t96 * t22) * s56 + t22 * t5 + t68 - t64 * (t20 * t3 -
     & (-(t96 + t82) * s26 - t108 + t8) * s26 - t32) * t36 - t1 * t100 *
     & t19 - t112 * t20 + t87 * (t70 - (-t106 - (t66 + t81) * s25 - t75)
     & * s26 - t7) * s34 - t104
      t82 = s26 + s56
      t101 = s25 * t44 - t58
      t102 = ddreal(10) * s16
      t105 = s16 * t24
      t106 = s16 * t13
      t108 = ddreal(9) * s25
      t112 = (-s16 + t108) * s25 + t8
      t114 = s16 - t45
      t117 = -t4 * t114
      t118 = t29 * s25
      t124 = s16 + t45
      t103 = s16 * t103
      t125 = ddreal(8) * t7
      t94 = t94 * t63
      t126 = t17 * t7
      t127 = (t66 - t83) * s25 + t46
      t128 = t17 * t21
      t129 = t39 * ((t102 - t61) * s25 + t5)
      t130 = t103 - t67
      t131 = (s25 * t130 + t110) * s25
      t132 = t39 * t57
      t133 = ddreal(18) * s25
      t83 = ((t3 * (-s34 * t13 * t63 + t7) + t100 * t36 - t24 - t35 * t1
     &30) * s56 + s25 * (((t120 - t61) * s25 - ddreal(22) * t5) * s25 + 
     &t125) + t36 * (s34 * t117 + (t96 + t42) * s16 + ddreal(8) * s25 * 
     &t97 - ddreal(11) * t21) + t6 + ddreal(4) * t132 * s26 - t4 * (t131
     & + t73 + t94) * s34) * s56 + t17 * (t132 * t18 + t98) + t3 * ((t6 
     &- t22) * s25 + t129 * s26) + (((-s34 * t124 + s16 * (-t43 + t56) +
     & (-s16 + t108 + t42) * s25) * s34 + (-ddreal(15) * t21 + t116) * s
     &25 + t17 * (t121 + t79) - t3 * (-s26 * t127 + t7)) * s34 + (((-t11
     &5 + t83) * s25 - t84) * s25 - ddreal(14) * t7) * s25 + (-t17 * (t1
     &31 + t73) - t119 * t13 * t63) * s26 + t6) * s34 + t105 * (-ddreal(
     &24) * s16 + t133)
      t84 = t64 * (-t56 * t24 + ((-s26 * t100 - ddreal(5) * t106 + ddrea
     &l(9) * t21) * s26 + ((-s16 + t81) * s25 - t15) * s25 + t99) * s26 
     &+ t32 * t13) * t36
      t97 = t87 * (-t3 * (t13 * t80 + t105) + t5 * ((-t26 + t45) * s25 +
     & t5) - (-ddreal(4) * s25 * t58 + ddreal(4) * t11 * t21 + ddreal(4)
     & * t7) * s26) * s34
      t83 = t61 * t39 * t101 * t18 + (-t37 * (s26 * t124 + t21) + t21 * 
     &t39 * ((-s25 + t26) * s25 + t8)) * s26 + (t83 * s56 + t3 * (t129 *
     & t18 + t21 * t6) + t4 * s34 * (t36 * (t24 * t4 + (-s26 * t114 * t3
     & + t112) * s26 + t35 * t13) + t64 * ((t21 * (-t103 + t45) - t125) 
     &* s25 + (-t3 * (((t96 - t61) * s25 + t75) * s25 + t73) - t94) * s2
     &6 + t6)) - t21 * t22 - t36 * (t36 * ((t28 + t11) * s26 + t21) - t2
     &1 * (t110 - t128) - ddreal(4) * t35 * t58 - s26 * (t18 * t113 + ((
     &s16 * t33 - ddreal(30) * s25) * s25 + ddreal(28) * t5) * s25 - t12
     &6) - t77 * t127) + t47 * s25 * t101 + ddreal(4) * t53 * t57 + s16 
     &* t22 * (-t102 + t109)) * s56 + t38 * (t118 * (t21 * t4 + t106) - 
     &t18 * (-s26 * t117 - t112) - t105) - t84 + t97 + t105 * ((-t5 - t2
     &) * s25 + t7) + t132 * t19 + t129 * t20
      t84 = s26 * t13
      t94 = (t107 + t45) * s25 - ddreal(51) * t5
      t96 = ddreal(15) * t6
      t97 = (s25 * t94 + t99) * s25 + t96
      t99 = -ddreal(49) * t5
      t100 = ddreal(21) * t7
      t101 = ddreal(5) * t6
      t56 = (((t81 - t56) * s25 + t99) * s25 + t100) * s25 + t101
      t102 = ddreal(10) * t5
      t106 = (-ddreal(21) * s16 + t81) * s25 + t102
      t108 = ddreal(14) * s16
      t110 = ddreal(26) * t5
      t67 = t26 + t67
      t112 = ddreal(15) * s25
      t113 = ddreal(15) * s26
      t93 = t93 * (s16 - t109)
      t117 = ddreal(47)
      t119 = (ddreal(192) * s16 + ddreal(81) * s25) * s25 + ddreal(78) *
     & t5
      t120 = ddreal(72) * s25
      t121 = ddreal(88) * s16
      t124 = t120 + t121
      t125 = ddreal(9) * s26
      t127 = t25 * t106
      t129 = ddreal(5) * t57 * t114
      t130 = ((ddreal(46) * s16 + ddreal(21) * s25) * s25 + ddreal(66) *
     & t5) * s25
      t131 = ddreal(19) * t7
      t132 = ddreal(23) * s16
      t134 = (t132 + t133) * s25 + t110
      t135 = ddreal(10) * t67
      t114 = t57 * t114
      t136 = ddreal(10) * t114 * s26
      t106 = t13 * t106
      t137 = ddreal(30) * s26
      t85 = s56 * (s56 * (s25 * (-t103 - t86 + t45) - t123 + t85) + s16 
     &* (s25 * (-t107 - t112) + t102) + t17 * (-t38 + t24) - s34 * (-s34
     & * (t113 + t115 + t133) + s16 * (ddreal(26) * s16 + t137) + s25 * 
     &(ddreal(25) * s26 + t132 + t133)) + t129 * s26) + s25 * (s25 * (t9
     &9 + t128) + t100) + s26 * (t127 + t136) - t3 * (t38 * (ddreal(5) *
     & t13 + t76 - s34) + t105) + t101 - s34 * (-s34 * (s25 * (ddreal(64
     &) * s16 + ddreal(27) * s25) + s26 * (t120 + t121 + t137) + t110) +
     & t130 + t131 + s26 * (s26 * t135 + ddreal(4) * t134))
      t85 = s56 * t85 + s25 * (t21 * t94 + t96) + t18 * (t106 * t17 + t1
     &36) + t3 * (s26 * t56 + t98) + ((((-s16 + t45 + t125) * s34 + (-t1
     &32 - t109) * s25 - ddreal(36) * t18 + t8 - ddreal(45) * t84) * s34
     & + s16 * ((ddreal(58) * s16 + ddreal(60) * s25) * s25 - t8) + (ddr
     &eal(30) * t18 + t119) * s26 + t33 * (t18 * t124 / ddreal(8) + t24)
     &) * s34 + (((-s16 * t117 - t109) * s25 - ddreal(58) * t5) * s25 - 
     &ddreal(52) * t7) * s25 + ((-s26 * t135 - t134 * t17) * s26 - t3 * 
     &(t130 + t131)) * s26 + t6) * s34
      t86 = s25 * t97 - t37 * (s16 - t95 - t45)
      t94 = t64 * (t19 * t3 - ((-(t115 + t112) * s26 - (s25 * t33 + ddre
     &al(42) * s16) * s25 - t110) * s26 - s16 * ((t133 + t71) * s25 - t8
     &)) * s26 - t122 + t51 * t4) * t36
      t6 = t87 * (t32 * (t108 + t28) - (-s16 * ((t108 + t133) * s25 + dd
     &real(19) * t5) + (-t67 * s26 - (t107 + t109) * s25 - t110) * s26) 
     &* s26 - t6) * s34
      t95 = t88 * (s25 * (t45 - t30) + t46)
      t98 = -s16 + s34
      t99 = ddreal(1) / t49
      t100 = ddreal(1) / s16
      t90 = ddreal(1) / (t90)**2
      t55 = ddreal(1) / (t55)**2
      t6 = s34 * (t95 + s26 * (-t38 * ((t113 * t13 - t93) * s26 + t17 * 
     &(t20 + t79) - t32 * t4) + t79 * t57 * ((-t81 - t44) * s25 + t46)) 
     &+ s56 * (s56 * t85 + t118 * t97 - ddreal(17) * t104 + t26 * t23 + 
     &t96 * t21 - t116 * t22 + t129 * t19 + t127 * t20 + t77 * t56 + (((
     &-t54 * (-t28 + t11 - t125) + t18 * (-ddreal(24) * s26 - ddreal(45)
     & * t13) + t4 * (s26 * t93 + t32) - t34) * s34 + ddreal(15) * t19 +
     & ((s26 * t124 + t119) * s26 + ((ddreal(120) * s16 + ddreal(24) * s
     &25) * s25 + ddreal(116) * t5) * s25 - t126) * s26 + t35 * ((ddreal
     &(34) * s16 + t133) * s25 - t75)) * s34 - t64 * (s16 * (((t115 + t1
     &33) * s25 + t117 * t5) * s25 - t73) + s26 * (((s25 * t117 + ddreal
     &(62) * s16) * s25 + ddreal(104) * t5) * s26 + ((ddreal(76) * s16 +
     & t111) * s25 + ddreal(94) * t5) * s25 + ddreal(57) * t7) + ddreal(
     &5) * t67 * t20)) * s34) + t18 * t86 + t20 * t56 + t94 + t114 * s26
     & * t19 + t106 * t19 - t6)
      t6 = -t9 * MBL1001(0) + (t6 * MBL1110(1) + t82 * t83 * MBL1011(1))
     & * t100 + t89 * t60 * t18 + (-t38 * (-t3 * t41 + t18 + t21 - t59) 
     &+ t21 * t1 * t13 * t57) * s26 + (t62 * s56 - s26 * (-t1 * (t60 * t
     &81 + ddreal(4) * t80) + t37 * t4) + t3 * (t18 * t78 - t58 * t79) -
     & ((((s25 - t11) * s25 - t17 * t41 + t77) * s34 + s25 * ((-t61 - t3
     &0) * s25 + t75) + ((s25 * t74 - t30 - t76) * s26 + t17 * (-ddreal(
     &4) * t14 + t5)) * s26) * s34 + t64 * (((t66 + t61) * s25 + ddreal(
     &15) * t5) * s25 + (ddreal(4) * t69 + (t71 + t72) * s25 + t15) * s2
     &6 - t73)) * s34 + t24 * t58) * s56 - t68 + t64 * (t18 * (t3 * (s16
     & + t61) + t29) + (s25 * t30 + t3 * t58) * s26 - t32) * t36 + t1 * 
     &t63 * t19 + t78 * t20 - t87 * (t70 - (-s25 * t66 - t3 * t65 - t69)
     & * s26 - t7) * s34 + t88 * t1
      t2 = t100 * (t55 * (-t82 * (-t18 * t98 - t4 * (s56 * (s26 * t98 - 
     &t35) - t35 * s26) - t52 * t98 + t79) * MBL1111D2(0) + s25 * s34 * 
     &(s26 * (-s26 * t91 + (-t4 * t84 - t35 + t54) * s34) + ((-t4 * t48 
     &- t2 + t36 + t5) * s56 + t4 * s26 * ((-s16 + s25) * s25 + t36 + t5
     &) - t35 * (-s16 - s25 + s34) - ddreal(4) * t54 * t13) * s56 + t51 
     &+ t35 * t13 * s26) * MBL1101(0) * t99) - s34 * (s25 + s26 + s56) *
     & MBL1110(0) * t90)
      result = -t55 * t99 * (-(t3 * (t1 * t51 - t50) - t4 * (s56 * t49 *
     & t52 + t53) + (t54 * ((t28 - t29 + t30) * s26 + t27) + t18 * (s26 
     &* t25 + t16) - t27 * t12 * s26 - t8 * t21) * s34 + ((s25 * ((t44 -
     & t45) * s25 - t46) + t3 * (-t38 + t7) + s34 * ((-t43 + t28 + t30) 
     &* s34 + t42 * t13 + t14 - t15) - t47) * s56 + (t4 * t16 * s26 - t3
     &3 * (-t31 + t32) - t34) * s34 - t17 * ((s26 * t39 + t38) * s26 - t
     &40) + t4 * (((t28 + t30) * s26 + t3 * (-t18 + t35)) * t36 + t41 * 
     &t10)) * s56 + t1 * t10 * t18 + t43 * t40) * t100 - t92 * I300s16(1
     &) + t6 * t90) / ddreal(4) + t2 / ddreal(2)

           intHLs16s25s26s34s56x1321D6eps0 = result
       end function intHLs16s25s26s34s56x1321D6eps0

       function intHLs16s25s26s34s56x1321D6eps1()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1321D6eps1
           type(dd_complex) ::  t1,t2

           type(dd_complex) :: result

      t1 = s25 + s26 - s34 + s56
      t2 = -ddreal(1) / ddreal(4)
      t1 = ddreal(1) / (t1)**2
      result = ddreal(1) / s16 * t1 * (t2 * (s25 + s26 + s34 + s56) - s3
     &4 * (s25 + s26 + s56) * MBL1110(1) / ddreal(2))

           intHLs16s25s26s34s56x1321D6eps1 = result
       end function intHLs16s25s26s34s56x1321D6eps1

       function intHLs16s25s26s34s56x1411D6eps0()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1411D6eps0
           type(dd_complex) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           type(dd_complex) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           type(dd_complex) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           type(dd_complex) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           type(dd_complex) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           type(dd_complex) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           type(dd_complex) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           type(dd_complex) ::  t175,t176,t177,t178,t179,t18,t180,t181,t182,t183,t184,t185
           type(dd_complex) ::  t186,t187,t188,t189,t19,t190,t191,t192,t193,t194,t195,t196
           type(dd_complex) ::  t197,t198,t199,t2,t20,t200,t201,t202,t203,t204,t205,t206
           type(dd_complex) ::  t207,t208,t209,t21,t210,t211,t212,t213,t214,t215,t216,t217
           type(dd_complex) ::  t218,t219,t22,t220,t221,t222,t223,t224,t225,t226,t227,t228
           type(dd_complex) ::  t229,t23,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239
           type(dd_complex) ::  t24,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249,t25
           type(dd_complex) ::  t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t26,t260
           type(dd_complex) ::  t261,t262,t263,t264,t265,t266,t267,t268,t269,t27,t270,t28
           type(dd_complex) ::  t29,t3,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39
           type(dd_complex) ::  t4,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t5
           type(dd_complex) ::  t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t6,t60
           type(dd_complex) ::  t61,t62,t63,t64,t65,t66,t67,t68,t69,t7,t70,t71
           type(dd_complex) ::  t72,t73,t74,t75,t76,t77,t78,t79,t8,t80,t81,t82
           type(dd_complex) ::  t83,t84,t85,t86,t87,t88,t89,t9,t90,t91,t92,t93
           type(dd_complex) ::  t94,t95,t96,t97,t98,t99

           type(dd_complex) :: result

      t1 = ddreal(3)
      t2 = t1 * s25
      t3 = -t2 + s16
      t4 = s16 - s25
      t5 = s16 + s25
      t6 = ddreal(6)
      t7 = ddreal(5)
      t8 = ddreal(4)
      t9 = (s16)**2
      t10 = (t9)**2
      t11 = t9 * t10
      t12 = s16 * t10
      t13 = s16 * t9
      t14 = t8 * s25
      t15 = t7 * s16
      t16 = t6 * t9
      t17 = (s25)**2
      t18 = (t17)**2
      t19 = s25 * t17
      t20 = t17 * t18
      t21 = s25 * t18
      t22 = -t21 + t12
      t23 = ddreal(33)
      t24 = ddreal(49) * s16
      t25 = t1 * t13
      t26 = s16 * s25
      t27 = t22 * t6 + t26 * (((-s25 * t23 - t24) * s25 - ddreal(41) * t
     &9) * s25 + t25)
      t28 = t18 + t10
      t29 = ddreal(42)
      t30 = ddreal(30)
      t31 = ddreal(2)
      t32 = t29 * t9
      t33 = t30 * t10
      t34 = (((ddreal(18) * t17 + t32) * s25 - ddreal(50) * t13) * s25 -
     & t33) * t17
      t35 = t31 * t11
      t36 = ddreal(9) * t26 * t28
      t37 = t1 * s16
      t38 = ddreal(18) * s25
      t39 = ddreal(53) * t9
      t40 = ddreal(9) * t13
      t41 = t7 * t10
      t42 = ddreal(16)
      t43 = s25 * t6
      t44 = t9 * t42
      t45 = ddreal(7) * t13
      t46 = ((-t43 + t37) * s25 - t44) * s25 + t45
      t47 = s25 + s26
      t48 = ddreal(37)
      t49 = ddreal(15)
      t50 = ddreal(12)
      t51 = s16 * t50
      t52 = t31 * t13
      t53 = -t19 + t13
      t54 = ddreal(7) * s16
      t55 = ddreal(9) * s25
      t56 = ddreal(18) * t13
      t57 = ddreal(28)
      t58 = ddreal(25)
      t59 = t6 * s16
      t60 = ddreal(19) * t9
      t61 = ((s25 * t58 - t59) * s25 - t60) * s25 + t13 * t57
      t62 = (s26)**2
      t63 = (t62)**2
      t64 = s26 * t63
      t65 = s26 * t62
      t66 = -t17 + t9
      t67 = t9 * t17
      t68 = t50 * t12
      t69 = ddreal(140) * s25
      t70 = ddreal(98) * s16
      t71 = t49 * t9
      t72 = t31 * t10
      t73 = ddreal(102) * s25
      t74 = ddreal(51) * t10
      t75 = ddreal(95) * s16
      t76 = ddreal(64) * t9
      t77 = ddreal(94) * t13
      t78 = ddreal(31)
      t79 = ddreal(39) * s25
      t80 = (s16 * t78 + t79) * s25 + t32
      t81 = t13 * t17
      t82 = t9 * s25
      t83 = t42 * t13
      t84 = ddreal(75) * s25
      t85 = ddreal(84) * t13
      t86 = (((ddreal(188) * s16 + t84) * s25 + ddreal(136) * t9) * s25 
     &+ t85) * s25 + t41
      t87 = ((ddreal(363) * s16 + ddreal(290) * s25) * s25 + ddreal(273)
     & * t9) * s25
      t88 = t87 + t77
      t89 = ddreal(84) * s25
      t90 = (ddreal(87) * s16 + t89) * s25 + ddreal(58) * t9
      t91 = ddreal(27) * s25
      t92 = t57 * s16
      t93 = -t91 - t92
      t94 = t26 * (((ddreal(152) * s16 + ddreal(160) * s25) * s25 + ddre
     &al(89) * t9) * s25 + t83)
      t95 = ddreal(36)
      t96 = s25 * t95
      t97 = t42 * s16
      t98 = ddreal(79) * s16
      t99 = ddreal(24) * t9
      t100 = ddreal(150) * t9
      t101 = ddreal(20) * t13
      t102 = ((ddreal(239) * s16 + ddreal(130) * s25) * s25 + t100) * s2
     &5
      t103 = t102 + t101
      t104 = (ddreal(55) * s16 + ddreal(50) * s25) * s25 + t99
      t105 = ddreal(69) * s16
      t106 = ddreal(82) * s25
      t107 = t106 + t105
      t108 = t26 * ((ddreal(100) * s25 + t98) * s25 + t99)
      t109 = ddreal(14) * s25
      t110 = t31 * s25
      t111 = t110 + s16
      t112 = ddreal(66) * s25
      t113 = ddreal(20) * t9
      t114 = (ddreal(78) * s16 + t112) * s25
      t115 = t114 + t113
      t116 = ddreal(51) * s16
      t117 = ddreal(73) * s25
      t118 = t117 + t116
      t119 = ddreal(10)
      t120 = t119 * s16
      t121 = t26 * t8
      t122 = ddreal(27) * s16
      t123 = ddreal(19) * t13
      t124 = ddreal(11) * t10
      t125 = t31 * t12
      t126 = ((((-t38 - t122) * s25 - ddreal(69) * t9) * s25 - t123) * s
     &25 + t124) * s25 + t125
      t127 = ddreal(187) * t13
      t128 = ddreal(54) * s16
      t129 = t18 * t5
      t130 = ddreal(39) * t9
      t131 = t13 * s25
      t132 = t7 * t65
      t133 = ddreal(51) * s25
      t134 = (((ddreal(837) * s16 + ddreal(705) * s25) * s25 + ddreal(46
     &8) * t9) * s25 + ddreal(561) * t13) * s25 + ddreal(153) * t10
      t135 = ddreal(63) * s16
      t136 = ddreal(47) * t13
      t137 = ((ddreal(86) * s25 + t135) * s25 + t39) * s25 + t136
      t138 = t12 * s25
      t139 = ddreal(8)
      t140 = ddreal(27) * t9
      t141 = t137 * t139
      t142 = t7 * s26
      t143 = t142 * t80
      t144 = s16 * t17
      t145 = t90 * t139
      t146 = t88 * t1
      t147 = ddreal(35) * s26
      t148 = t104 * t50
      t149 = ddreal(64) * s26
      t150 = ddreal(32) * s25
      t151 = t23 * s26
      t128 = ((((-((ddreal(20) * s16 + t96 + t151) * s26 + t121) * s34 +
     & ((t1 * t118 + t149) * s26 + t8 * t115 / ddreal(2)) * s26 + t26 * 
     &(t97 + t150)) * s34 - (((t107 * t8 + t147) * s26 + t148) * s26 + t
     &31 * t103) * s26 - t108) * s34 + (((-t142 * t93 + t145) * s26 + t1
     &46) * s26 + t8 * t86) * s26 + t94) * s34 - t62 * ((t143 + t141) * 
     &s26 + t134) - t8 * ((((((ddreal(162) * s16 + t133) * s25 + ddreal(
     &95) * t9) * s25 + ddreal(54) * t13) * s25 + ddreal(57) * t10) * s2
     &5 + t12) * s26 + t138) - t144 * (((ddreal(142) * s16 + t69) * s25 
     &+ t140) * s25 + ddreal(75) * t13)) * s34 + t47 * (ddreal(64) * t12
     &9 * s16 + s26 * (((((ddreal(227) * s25 + t128) * s25 - t9) * s25 -
     & ddreal(68) * t13) * s25 + ddreal(144) * t10) * s26 + ((((ddreal(2
     &28) * s16 + ddreal(76) * s25) * s25 + t76) * s25 - t127) * s25 + d
     &dreal(66) * t10) * s25 + ddreal(45) * t12) + t131 * ((-s25 * t48 -
     & t128) * s25 + t130) + t132 * t61)
      t152 = (s34)**2
      t153 = (t152)**2
      t154 = s34 * t152
      t155 = t154 * t153
      t156 = t152 * t153
      t157 = s34 * t153
      t39 = t4 * ((((-t38 - t37) * s25 - t39) * s25 + t40) * s25 + t41)
      t158 = t39 * t62
      t159 = t50 * s25
      t160 = ddreal(11) * s16
      t161 = ddreal(17) * t13
      t162 = ddreal(13) * t10
      t163 = t12 * t6
      t164 = t1 * t4
      t165 = ddreal(21)
      t166 = t165 * t9
      t167 = (ddreal(94) * s16 + ddreal(101) * s25) * s25
      t168 = ddreal(85) * t13
      t169 = ddreal(70) * t10
      t170 = t49 * t12
      t171 = ddreal(44) * s25
      t172 = t50 * t9
      t173 = ddreal(18) * t10
      t174 = (((t171 + t37) * s25 - t172) * s25 + t40) * s25 + t173
      t175 = ddreal(112) * t13
      t176 = t10 * t6
      t177 = t12 * t29
      t178 = s26 * t119
      t179 = t178 * t61
      t180 = ddreal(108) * t13
      t181 = t178 * t80
      t182 = -t119 * t93 * s26
      t183 = ddreal(70) * s26
      t184 = t12 * t17
      t185 = t139 * t9
      t186 = s16 * t19
      t187 = t186 * (((t97 + t91) * s25 - t185) * s25 - ddreal(44) * t13
     &)
      t188 = t119 * t4 * t46
      t189 = -ddreal(88)
      t190 = ddreal(168) * s25
      t191 = ddreal(116) * s16
      t192 = t4 * t46
      t193 = s16 * t49
      t194 = t153 + t10
      t195 = ddreal(23)
      t196 = ddreal(9) * s16
      t197 = (-t196 + t43) * s25
      t198 = t195 * t13
      t199 = (t47)**2
      t105 = (s25 * ((t197 + t60) * s25 - t198) + ddreal(7) * t194 + ((s
     &34 * t93 + t80) * s34 - t61) * s34) * s56 + s25 * ((((-t193 + t38)
     & * s25 + ddreal(50) * t9) * s25 - ddreal(62) * t13) * s25 + t10 * 
     &t8) - t157 * t42 + t7 * (s26 * t192 + t12) + ((((t147 + t106 + t10
     &5) * s34 - s16 * (ddreal(140) * s26 + t191) - (ddreal(174) * s16 +
     & ddreal(135) * s26 + t190) * s25) * s34 + ((ddreal(126) * s16 + dd
     &real(172) * s25) * s25 + ddreal(106) * t9) * s25 + t143 + t77) * s
     &34 + (((s25 * t189 - t59) * s25 + t99) * s25 - t56) * s25 - t10 * 
     &t95 - t142 * t61) * s34
      t87 = t105 * s56 + t34 + t35 + t36 + ((((-(t117 + t116 + t149) * s
     &34 + (ddreal(220) * s16 + ddreal(200) * s25) * s25 + (ddreal(276) 
     &* s16 + ddreal(328) * s25 + t183) * s26 + ddreal(96) * t9) * s34 -
     & (t145 + t182) * s26 - t77 - t87) * s34 + (((ddreal(279) * s16 + d
     &dreal(235) * s25) * s25 + ddreal(156) * t9) * s25 + t127) * s25 + 
     &(t181 + t141) * s26 + t74) * s34 + (((-t167 - t166) * s25 + t168) 
     &* s25 - t169) * s25 + (-t139 * t174 - t179) * s26 - t170) * s34 + 
     &ddreal(11) * t156 + t39 * t8 * s26 + t188 * t62
      t87 = t87 * s56 + t1 * (s26 * t126 * t4 - t184) - t155 * t31 + t6 
     &* ((t20 + t11) * s25 + t158) + ((((((t38 + t120 + t151) * s34 - (d
     &dreal(153) * s16 + ddreal(219) * s25 + ddreal(96) * s26) * s26 - t
     &113 - t114) * s34 + ((t107 * t6 + t183) * s26 + t148) * s26 + t101
     & + t102) * s34 - (((ddreal(376) * s16 + ddreal(150) * s25) * s25 +
     & ddreal(272) * t9) * s25 + ddreal(168) * t13) * s25 - ((t50 * t90 
     &+ t182) * s26 + t146) * s26 - t10 * t119) * s34 + ((((ddreal(324) 
     &* s16 + t73) * s25 + ddreal(190) * t9) * s25 + t180) * s25 + ddrea
     &l(114) * t10) * s25 + ((t137 * t50 + t181) * s26 + t134) * s26 + t
     &125) * s34 + s25 * (((((-ddreal(146) * s16 - ddreal(38) * s25) * s
     &25 - t76) * s25 + t175) * s25 - t176) * s25 - t177) + ((-t174 * t5
     &0 - t179) * s26 - t1 * ((((t167 + t166) * s25 - t168) * s25 + t169
     &) * s25 + t170)) * s26) * s34 + t187 + t188 * t65
      t69 = t47 * (t50 * t81 * t5 + (s16 * ((((t69 + t70) * s25 + t71) *
     & s25 + ddreal(63) * t13) * s25 + t72) + ((t80 * s26 + ((ddreal(133
     &) * s25 + t75) * s25 + t76) * s25 + t77) * s26 + (((ddreal(184) * 
     &s16 + t73) * s25 + ddreal(92) * t9) * s25 + ddreal(93) * t13) * s2
     &5 + t74) * s26) * s26 + t82 * (ddreal(44) * t19 + t52)) * t152
      t56 = t199 * (-ddreal(26) * t67 * t66 + s26 * (s16 * ((((ddreal(64
     &) * s25 + t51) * s25 - t48 * t9) * s25 - t52) * s25 + t10 * t49) +
     & t61 * t62) + t31 * (s16 * (t17 * (t55 + t54) + t56) - ddreal(19) 
     &* s25 * t53) * t62 + t68 * s25) * s34
      t61 = t4 * s25
      t73 = t13 * t19
      t74 = t155 * t62 - t53 * t73
      t77 = t9 * t18
      t27 = -t31 * t74 - t6 * t77 * t53 + (t156 * ((ddreal(11) * s26 + t
     &38 + t120) * s26 + t121) + t144 * t3 * t4 * t5 * ((t15 + t14) * s2
     &5 + t16)) * s26 + (t87 * s56 - s34 * t128 + t13 * t21 + t8 * s26 *
     & (t158 - t155) + t144 * (((t17 * (t160 + t159) - t161) * s25 - t16
     &2) * s25 + t163) + t110 * t4 * t27 * s26 + t164 * t126 * t62 + t7 
     &* t4 * t46 * t63) * s56 + (t36 + t35 + t34) * t65 + ((-(t42 * s26 
     &* (t111 * t26 + t65) + t62 * (s26 * t118 + t115) + t67 * t31) * s3
     &4 + t65 * (t104 * t8 + ddreal(7) * t62) + ((t107 * t62 + t103) * s
     &26 + t108) * s26 + t67 * (t109 + t59)) * s34 - t31 * t62 * (t62 * 
     &t90 + t86) + (t62 * (t62 * t93 - t88) - t94) * s26 + t67 * ((-t97 
     &- t96) * s25 - t16)) * t154 + t69 - t56 + t192 * t64 + t39 * t63 +
     & t61 * t27 * t62
      t34 = -t31 * (s34 * t5 + t26) + t152 + t17 + t9
      t35 = s25 + s26 - s34 + s56
      t36 = t31 * s26
      t39 = (-s34 + t5 + s26) * s26 + (t36 + s16 + s25 - s34 + s56) * s5
     &6 + t26
      t46 = t7 * t9
      t56 = t10 * t5
      t69 = ddreal(14) * s16
      t74 = t1 * t129
      t80 = t67 * (t43 - t69)
      t86 = t82 * t31
      t53 = -t86 + t53
      t87 = t119 * s25
      t88 = t87 + t59
      t90 = t8 * t13
      t93 = ddreal(13) * s25
      t94 = ddreal(18) * s16
      t96 = t119 * t13
      t102 = (t17 * (t94 + t93) + t96) * s25 + t41
      t103 = t19 + t13
      t104 = t31 * t103 - t82
      t105 = t6 * s26
      t106 = t1 * t10
      t107 = ddreal(17) * s16
      t108 = ddreal(24) * s25
      t114 = ddreal(32) * s16
      t115 = ddreal(22) * s25
      t117 = ((t115 + t114) * s25 + t99) * s25
      t118 = t117 + t96
      t125 = t17 + t9
      t126 = t125 * t6 + t26 * t7
      t127 = t26 * ((t108 + t107) * s25 + t172)
      t128 = t31 * s16
      t129 = t2 + t128
      t134 = (t55 + t160) * s25
      t137 = t134 + t46
      t141 = t50 * s26
      t143 = t1 * s26
      t145 = ddreal(7) * s25
      t146 = t26 * t31
      t148 = ddreal(9) * s26
      t149 = t148 * t126
      t151 = t5 * s26
      t158 = t31 * t154
      t160 = (s56)**2
      t167 = (t160)**2
      t168 = s56 * t167
      t169 = s56 * t160
      t174 = t9 * t21
      t109 = (t80 + t74 - t158 * (ddreal(18) * t151 + t46 + t134) + ((t1
     &52 * (t145 + t148 + t15 - s34) + t117 + t149 + t96) * s34 - (ddrea
     &l(13) * t19 + t96) * s25 - ddreal(18) * s26 * t104 - ddreal(18) * 
     &t186 - t41) * s34 + t56 + t148 * t4 * t53) * s56 + t31 * ((-s26 * 
     &t153 - s26 * t102 + t81 - t26 * (t17 * t88 + t90) - ddreal(9) * t1
     &04 * t62) * s34 + t138 + t4 * ((((-t59 - t2) * s25 - t172) * s25 +
     & t52) * s25 + t10) * s26) - t8 * (t154 * ((t148 * t5 + t137) * s26
     & + t26 * t129) + t73) + t152 * (t152 * ((t109 + t120 + t148) * s26
     & + t146) + (t118 * t31 + t149) * s26 + t127) + t144 * (t17 * (t43 
     &+ t37) - t45) + ddreal(9) * t4 * t53 * t62
      t117 = t26 * t4
      t134 = t13 * t18
      t149 = ddreal(22) * s16
      t179 = t8 * s16
      t181 = t1 * t28 + t26 * ((s25 * t195 - t179) * s25 - t46)
      t182 = t7 * t13
      t135 = (((t135 + t91) * s25 - ddreal(38) * t9) * s25 + t182) * s25
     & + t106
      t183 = ddreal(17) * t9
      t187 = ((t193 + t91) * s25 - t183) * s25 + t182
      t188 = t57 * s25
      t189 = ddreal(32) * t9
      t192 = (t59 + t188) * s25
      t200 = t111 * t65
      t201 = s25 * (t10 - t200)
      t202 = t165 * s25
      t203 = t29 * t13
      t204 = t49 * s25
      t205 = t30 * t9
      t206 = ddreal(13) * s16
      t207 = ddreal(89) * s16
      t208 = ((s25 * t42 + t207) * s25 + ddreal(90) * t9) * s25
      t209 = t208 + t96
      t210 = s16 * t48
      t211 = (t210 + t204) * s25
      t212 = t211 + t99
      t213 = t15 + t2
      t92 = t26 * ((t204 + t92) * s25 + t16)
      t214 = t213 * t8
      t215 = ddreal(11) * s25
      t216 = t119 * t9
      t116 = (ddreal(29) * s25 + t116) * s25
      t217 = t116 + t216
      t218 = t78 * s25
      t219 = t218 + t122
      t220 = ddreal(27) * t10
      t221 = s16 * t165
      t222 = ddreal(11) * t9
      t223 = ddreal(14) * t13
      t224 = ((-t159 + t37) * s25 + t60) * s25 + t223
      t225 = (-t55 + s16) * s25 + t71
      t226 = t10 * s25
      t227 = ddreal(34)
      t228 = t149 * t18
      t229 = t9 * t227
      t230 = s16 * t139
      t231 = ddreal(47) * t9
      t232 = ddreal(9) * t10
      t233 = t8 * t4 * t187
      t234 = t13 * t139
      t235 = ((ddreal(108) * s16 + ddreal(312) * s25) * s25 + ddreal(192
     &) * t9) * s25 - ddreal(120) * t13
      t236 = t9 * t19
      t237 = t29 * t10
      t238 = ddreal(60) * s25
      t239 = t30 * s25
      t240 = ddreal(81) * s16
      t241 = t30 * s26
      t242 = ddreal(93) * s25
      t243 = ddreal(18) * s26
      t244 = ddreal(24) * t13
      t245 = t3 * t65
      t246 = s25 * s26
      t247 = t246 + t9
      t248 = ddreal(40) * s26
      t249 = ddreal(20) * s26
      t250 = s16 * t95
      t251 = t6 * t157
      t252 = ddreal(60) * s26
      t253 = s26 * t49
      t91 = (s25 * (((t51 - t91) * s25 + t189) * s25 - ddreal(22) * t13 
     &- t253 * t3 * t4) + t194 * t7 + s34 * (((-t253 - t38 + t128) * s25
     & + t205) * s34 + (s16 * (t241 + t114) + (ddreal(52) * s25 + t94 + 
     &t252) * s25) * s25 - t101) - t214 * t154) * s56 + t1 * t12 + t17 *
     & (((-t250 - t91) * s25 + ddreal(101) * t9) * s25 - ddreal(43) * t1
     &3) + t31 * (-t154 * (s16 * t248 + t211 + ddreal(24) * t247) + t226
     &) + ((t139 * t225 * s26 + t152 * (t249 + t218 + t122) - t19 * t95 
     &- t239 * t62 + ddreal(9) * t144 + t203 + ddreal(57) * t82) * s34 +
     & ddreal(2) / ddreal(3) * t235 * s26 - t173 + ddreal(68) * t18 + dd
     &real(40) * t67 - t131 * t42 + ddreal(74) * t186 + t238 * t111 * t6
     &2) * s34 + t233 * s26 - t61 * t30 * t3 * t62 - t251
      t22 = t91 * s56 + ddreal(9) * s25 * t22 + t1 * (s26 * t135 * t4 + 
     &t73) + ((t50 * t225 * t62 + (((-t243 - t159 - t15 + s34) * s34 + (
     &t242 + t241 + t240) * s26 + t116 + t216) * s34 - (ddreal(24) * s26
     & * t213 + t212 * t6) * s26 - t208 - t96) * s34 + t148 * t224 - t16
     &5 * t18 - t186 * t49 - t239 * t65 + ddreal(84) * t131 + t41 + ddre
     &al(57) * t67) * s34 + ddreal(118) * s16 * t18 - t237 * s25 + t200 
     &* t238 + t21 * t57 + t235 * t62 + t236 * t31 - t12 + ddreal(27) * 
     &t81 - t105 * ((((-s25 * t227 - t210) * s25 - t113) * s25 + t234) *
     & s25 + t232)) * s34 + t144 * (t17 * (t240 - t238) - t244) - t245 *
     & t61 * t30 + t6 * t4 * t187 * t62
      t22 = t22 * s56 + t31 * ((t152 * (t152 * (s26 * s34 - (t159 + t148
     & + t15) * s26 - t26) - t92 - t139 * t213 * t65 - s26 * (s26 * t1 *
     & t212 + t209)) + t47 * (-t201 * t49 - (-(((ddreal(74) * s25 + t221
     &) * s25 + t76) * s25 - ddreal(40) * t13) * s26 - (((ddreal(90) * s
     &16 + t188) * s25 - t8 * t9) * s25 + t83) * s25 + t220) * s26 - t12
     & + t144 * (t192 + t222))) * s34 + t174) + t152 * (t152 * (ddreal(2
     &0) * t65 + s26 * (s26 * t1 * t219 + t217 * t31) + t26 * (t230 + t1
     &15)) + t139 * (t225 * t65 + t226) - t228 + ddreal(9) * t224 * t62 
     &- t204 * t63 + t36 * ((((-t202 - t193) * s25 + ddreal(57) * t9) * 
     &s25 + t85) * s25 + t41) + t229 * t19 + ddreal(55) * t81) + t144 * 
     &(((-ddreal(24) * t17 + t231) * s25 - t13 * t227) * s25 + t232) + t
     &105 * t61 * t181 - t61 * t49 * t3 * t63 + t164 * t135 * t62 + t233
     & * t65
      t22 = t22 * s56 + t153 * ((t132 + t146 * (t215 + t179)) * s26 + t6
     &2 * (s26 * t219 + t217) + t67) + t154 * (-t31 * s26 * (t212 * t62 
     &+ t92) + t62 * (-t214 * t62 - t209) - t67 * (t87 + t37)) + t156 * 
     &t62 + (t62 * (s26 * t187 + t135) + t144 * (((t108 + t149) * s25 - 
     &t58 * t9) * s25 + t40)) * s26 * t4 + t199 * ((ddreal(18) * t186 - 
     &t36 * (-t117 * t50 - t19 * t57 + t40)) * s16 - t201 * t6 - t125 * 
     &t13 + t62 * ((t192 + t189) * s25 - t101)) * s34 + t47 * (-t1 * s25
     & * (t63 - t10) + (s16 * (((s16 * t30 - t115) * s25 + t32) * s25 + 
     &t182) + (((-t204 + t128) * s25 + t205) * s26 + ((-t202 + t54) * s2
     &5 + t140) * s25 + t203) * s26) * s26 + t67 * (t14 + t206)) * t152 
     &- s26 * ((t105 + t159 + t15) * s26 + t146) * t157
      t61 = s26 + s56
      t76 = ddreal(18) * t9
      t83 = t6 * t13
      t85 = (((t114 + t2) * s25 - t76) * s25 - t83) * s25 + t232
      t91 = -t30 * t144 * t4 + (ddreal(9) * t19 + t234) * s25 + t106
      t92 = ((t55 + t230) * s25 - t172) * s25 + t182
      t94 = t13 * t50
      t108 = t9 * t57
      t116 = t8 * t17
      t132 = t37 - t110
      t135 = t30 * t13
      t164 = t65 * t132
      t173 = t36 * s16
      t187 = ddreal(24) * s16
      t188 = ddreal(14) * t9
      t192 = t145 + t37
      t201 = (t204 + t250) * s25
      t208 = ((t201 + t231) * s25 + ddreal(45) * t13) * s25 + t41
      t209 = t103 * t6 + t26 * (t120 + t145)
      t211 = (t187 + t202) * s25 + t113
      t212 = t1 * t9
      t213 = t211 * s26
      t214 = ddreal(51) * t9
      t201 = (t201 + t214) * s25
      t217 = t201 + t96
      t219 = ddreal(17) * s25
      t224 = (t219 + t210) * s25 + t140
      t225 = t15 + t14
      t185 = t26 * ((t149 + t215) * s25 + t185)
      t226 = t37 + s25
      t231 = t226 * t62
      t233 = ((t51 + t2) * s25 + t46) * s26
      t234 = ddreal(123) * t9
      t235 = t29 * s25
      t240 = ddreal(105) * s16
      t254 = ddreal(58) * t13
      t255 = ddreal(84) * t10
      t256 = t7 * t12
      t257 = ((t79 + s16) * s25 + t188) * s25 + t135
      t258 = t138 * t139
      t259 = t8 * t257
      t99 = t1 * ((((ddreal(62) * s25 + t210) * s25 + t99) * s25 + ddrea
     &l(59) * t13) * s25 + t237)
      t260 = ddreal(55) * t10
      t261 = t7 * s25
      t262 = t261 * t63
      t263 = t1 * t224
      t264 = (t4)**2
      t265 = s16 * s34
      t266 = t246 * t264
      t267 = t152 * (t152 * (t185 + t8 * t225 * t65 + s26 * (s26 * t263 
     &+ t217 * t31)) + t259 * t65 + t99 * t62 + t258 + t36 * (((((t240 +
     & t239) * s25 - t188) * s25 + t254) * s25 + t255) * s25 + t256) - t
     &77 * t6 + ddreal(27) * t73 + t260 * t17 + ddreal(66) * s16 * t21 +
     & t261 * t192 * t63)
      t268 = ddreal(9) * t12
      t269 = t264 * t62
      t270 = t8 * t264
      t97 = t31 * ((t152 * (t152 * (s26 * t265 - t1 * t231 - t233 - t82)
     & - t228 - t176 * s25 - t81 * t57 - ddreal(17) * t236 - t50 * t209 
     &* t62 - t262 - s26 * (t211 * t31 * t62 + t208 * t31)) - t47 * (s25
     & * ((((-t76 + ddreal(22) * t26) * s25 - t123) * s25 - t124 - t164 
     &* t119) * s25 + t170) + (((((-ddreal(26) * s16 + t235) * s25 - t20
     &5) * s25 - ddreal(56) * t13) * s25 + ddreal(40) * t10) * s26 + (((
     &(ddreal(62) * s16 + t204) * s25 - t234) * s25 + t182) * s25 - t10 
     &* t42) * s25 + ddreal(27) * t12) * s26 + t11)) * s34 + t266 * t85)
     & + t267 + t144 * (((((-t97 + t215) * s25 - t183) * s25 + t136) * s
     &25 - t10 * t227) * s25 + t268) - t262 * t3 * t264 + t269 * t1 * t9
     &1 + t270 * t92 * t65
      t123 = ddreal(19) * s25
      t124 = t10 * t139
      t136 = t26 * t125
      t57 = t57 * t136
      t170 = ddreal(20) * t10
      t71 = t17 * (t17 * t78 - t71) - t57 + t170
      t78 = ddreal(84) * s16
      t176 = t87 * t192
      t182 = t264 * s26
      t183 = ddreal(24) * t10
      t205 = t269 * t6
      t28 = t1 * (t182 * t91 + t19 * t28) + ((t6 * t257 * t62 + t255 * s
     &25 + t99 * s26 + t152 * ((-s16 * (t243 + t15) - (t51 + t105 + t2) 
     &* s25 + t265) * s34 + (t105 * t225 + t263) * s26 + t201 + t96) + t
     &176 * t65 + t240 * t18 + t21 * t30 - ddreal(14) * t236 + t256 + dd
     &real(58) * t81) * s34 + (((((-t204 - t78) * s25 + ddreal(141) * t9
     &) * s25 + t223) * s25 + t220 + ddreal(40) * t164) * s25 - t177) * 
     &s25 - t6 * s26 * (t71 * s26 + ((((t123 + t51) * s25 - t214) * s25 
     &- t161) * s25 + t124) * s25 + t268) - t11) * s34 + t26 * ((t17 * (
     &(ddreal(26) * s25 - t98) * s25 + ddreal(62) * t9) - t183) * s25 + 
     &t268) - t158 * (t1 * t211 * t62 + t141 * t209 + t18 * t49 + t186 *
     & t95 + t65 * t87 + ddreal(45) * t131 + t41 + ddreal(47) * t67) - t
     &87 * t245 * t264 + t205 * t92
      t87 = ddreal(24) * t67
      t57 = s25 * ((((t55 - t120) * s25 - t60) * s25 + t13 * t48) * s25 
     &- ddreal(22) * t10) + t7 * (-t266 * t3 + t12) + (((s34 * t225 - (t
     &187 + t202 + t178) * s25 - t113) * s34 + (s16 * (t253 + t69) + (t7
     &9 + t147 + s16) * s25) * s25 + t135) * s34 + t17 * (s16 * (t193 + 
     &t252) + (-t218 - t248) * s25) - t170 + t57) * s34
      t38 = (s34 * (t62 * ((s26 * t225 + t224) * s26 + t217) + t81 + t18
     &5 * s26) - t31 * s26 * ((s25 * t65 + t208) * s26 + t26 * (((t115 +
     & t107) * s25 + t108) * s25 + t83)) + t65 * (-t139 * t209 - t213) +
     & t67 * ((-t120 - t43) * s25 - t212)) * s34 + t47 * (t82 * (((t38 -
     & t128) * s25 + ddreal(13) * t9) * s25 + t25) + (s16 * ((((t112 - t
     &187) * s25 + ddreal(29) * t9) * s25 + t203) * s25 + t41) + ((t246 
     &* t192 + ((t150 - t128) * s25 + t188) * s25 + t135) * s26 + (((ddr
     &eal(39) * s16 + t239) * s25 + t216) * s25 + ddreal(29) * t13) * s2
     &5 + t237) * s26) * s26)
      t38 = t38 * s34 - t199 * ((((-t135 + ddreal(18) * t82) * s25 - t16
     &4 * t8) * s25 + t163) * s25 + t10 * t125 + t62 * (-t49 * t17 * t66
     & + s16 * ((-t108 - t116) * s25 + t101)) + t173 * ((((-t250 + t115)
     & * s25 + t222) * s25 - t94) * s25 + t232))
      t60 = t144 * t264
      t69 = t264 * s25
      t79 = s16 * t156
      t95 = s25 * (s34 * t132 * t14 + t1 * t19 + t152 * t192 + t7 * t82 
     &- t13 - ddreal(7) * t144 - t158) * t168
      t98 = ddreal(19) * s16
      t99 = t236 * ((((-t98 + t43) * s25 + t195 * t9) * s25 - ddreal(13)
     & * t13) * s25 + t106)
      t101 = t128 + s25
      t107 = -t17 * t4 * t8 + t3 * t9
      t113 = (t43 + t179) * s25 + t9
      t115 = (t43 + s16) * s25 + t212
      t147 = s25 * t113
      t150 = t14 + t37
      t158 = t1 * t107
      t161 = t1 * t115
      t163 = ((-t36 * t144 * (-t2 + t128) + t62 * (s26 * t115 + t147) + 
     &t236) * s34 - t47 * (t67 * (-t179 + t110) + s26 * (s26 * t107 - t2
     &6 * (-t197 + t9)))) * s34 + t236 * t125
      t164 = t17 * t264
      t107 = s25 * (-t186 * t31 + (((-t150 + s34) * s34 + t115) * s34 - 
     &t107) * s34 + t125 * t17) * t169
      t107 = s25 * t163 - (t1 * t19 * (t67 - t269) - t31 * t19 * (t101 *
     & t182 + t186) - ((s25 * (s26 * t110 * t113 + t161 * t62 + t186 * t
     &6 - t67 * t8) + (-s25 * ((s26 * t1 * t150 + t111 * t14) * s26 + t1
     &44 * t31) + t246 * (t143 + t110) * s34) * s34) * s34 - s25 * (-s26
     & * t101 * t116 * t132 + t158 * t62 + t18 * t59 - t7 * t81 - ddreal
     &(7) * t236)) * s34 - t10 * t19) * s56 - t160 * (t1 * t19 * (t82 - 
     &t182) - s25 * (((t161 * s26 + ((t143 + s25) * s34 - (t141 + t14 + 
     &t128) * s25 - t148 * s16) * s34 + t147) * s34 - t17 * ((t14 + t128
     &) * s25 - t172) - t158 * s26) * s34 + t21) - t73 * t31) + s25 * (t
     &101 * t164 + t153 * t47) * t62 + t246 * (-t154 * (t150 * t62 + t31
     & * s25 * (s26 * t111 + t26)) + t60 * t111) + t107 + t264 * t19 * t
     &65
      t113 = s25 + s26 + s56
      t115 = -t261 + t37
      t16 = (t230 + s25) * s25 + t16
      t116 = t1 * t26 + t125
      t147 = (t230 + t43) * s25 + t9
      t148 = t1 * t4 * t9 + t19
      t150 = t148 * t65
      t158 = -ddreal(48)
      t161 = ddreal(117) * s25
      t163 = ddreal(20) * s25
      t170 = ddreal(60) * t13
      t172 = (t196 + t215) * s25 + t76
      t177 = t172 * s26
      t179 = ddreal(177) * s25
      t185 = ddreal(52) * t13
      t188 = ddreal(40) * s25
      t192 = ddreal(60) * t9
      t193 = t59 + t261
      t197 = t31 * t193
      t201 = ddreal(108) * t9
      t203 = ((ddreal(165) * s16 + t188) * s25 + t229) * s25
      t204 = t203 - t254
      t112 = (ddreal(106) * s16 + t112) * s25
      t208 = t112 + t166
      t209 = t10 * t17
      t211 = ddreal(11) * t73
      t212 = (ddreal(58) * s16 + t2) * s25 + t216
      t214 = (-t163 - t128) * s25
      t216 = t214 + t189
      t217 = s16 * t212
      t93 = t93 + t15
      t220 = -t54 + t14
      t223 = s16 * t93
      t224 = -ddreal(362)
      t100 = (((s16 * t224 - t179) * s25 - t100) * s25 + ddreal(106) * t
     &13) * s25 + t41
      t224 = -ddreal(393)
      t89 = ((ddreal(293) * s16 + t89) * s25 + ddreal(226) * t9) * s25
      t225 = t89 + t244
      t207 = (ddreal(54) * s25 + t207) * s25
      t227 = t207 + t192
      t228 = t227 * t7
      t230 = t8 * t225
      t232 = ((((s16 * t224 - t188) * s25 - ddreal(378) * t9) * s25 - dd
     &real(121) * t13) * s25 + ddreal(52) * t10) * t1
      t240 = (-t98 - t242) * s25
      t246 = t1 * t204
      t248 = t8 * t208
      t250 = t1 * t216
      t256 = s16 * s26
      t257 = t269 * t115
      t262 = t115 * t264
      t263 = t111 * t63
      t72 = (-t50 * t193 * t64 + ((s34 * (-t31 * s16 * (s26 * t93 + t26)
     & + t62 * (t1 * t220 + t141)) + t131 * t139 + t200 * t158 + t173 * 
     &t212 + t250 * t62 - t30 * t63 + t87) * s34 - t50 * t82 * t66 + (t6
     &5 * (ddreal(150) * t5 + t243) - t128 * ((t240 + t201) * s25 + t96)
     &) * s26 + t62 * (s26 * t248 + t246) - ddreal(65) * t81) * s34 + t1
     &73 * t100 - t228 * t63 - t230 * t65 + t232 * t62 + ddreal(70) * t2
     &09 + t258 - ddreal(69) * t73 - ddreal(150) * t77) * s34 + t47 * (t
     &9 * ((((ddreal(131) * s16 + t190) * s25 + ddreal(142) * t9) * s25 
     &- ddreal(41) * t13) * s25 - t72) + (((((ddreal(436) * s16 + ddreal
     &(144) * s25) * s25 + ddreal(282) * t9) * s25 + ddreal(300) * t13) 
     &* s26 + (((ddreal(708) * s16 + t238) * s25 + ddreal(514) * t9) * s
     &25 + ddreal(672) * t13) * s25 + t255) * s26 + s16 * ((((ddreal(698
     &) * s16 + ddreal(234) * s25) * s25 + ddreal(303) * t9) * s25 + ddr
     &eal(336) * t13) * s25 - ddreal(69) * t10)) * s26 + t6 * t172 * t63
     &)
      t73 = ddreal(38) * s16
      t93 = -t73 + t239
      t212 = t270 * t115
      t239 = t212 * t147
      t255 = ddreal(18) * t12
      t258 = t252 * t148
      t266 = t165 * t10
      t133 = (((ddreal(286) * s16 + t133) * s25 + ddreal(199) * t9) * s2
     &5 + ddreal(243) * t13) * s25
      t267 = t133 + t266
      t24 = ((t202 + t24) * s25 + t130) * s25 + t135
      t130 = t253 * t172
      t241 = t241 * t193
      t253 = ddreal(300) * t5
      t270 = ddreal(72)
      t100 = (((((-t221 + t243 + t159) * s26 - t223 + t265) * s34 + ((-t
     &111 * t270 - t252) * s26 + t250) * s26 + t217) * s34 + s16 * ((-t2
     &40 - t201) * s25 - t96) + (((ddreal(45) * s26 + t253) * s26 + t208
     & * t6) * s26 + t246) * s26) * s34 + s16 * t100 + (((-t119 * t227 -
     & t241) * s26 - t225 * t6) * s26 + t232) * s26) * s34 + t62 * ((ddr
     &eal(20) * t24 + t130) * s26 + t267 * t6) + (((((ddreal(942) * s16 
     &+ t238) * s25 + ddreal(1212) * t9) * s25 + ddreal(975) * t13) * s2
     &5 + ddreal(420) * t10) * s25 - ddreal(69) * t12) * s26 - t11 + t26
     & * ((((ddreal(433) * s16 + t161) * s25 + ddreal(217) * t9) * s25 +
     & ddreal(239) * t13) * s25 - t260)
      t221 = ddreal(57) * s25
      t225 = ddreal(109) * t9
      t187 = (-t187 + t221) * s25
      t12 = t50 * (t21 + t12)
      t232 = (((t98 + t159) * s25 - ddreal(46) * t9) * s25 + t40) * s25
      t238 = t232 + t33
      t240 = t182 * t115
      t206 = s25 * (-t261 + t206) - t222
      t32 = t199 * (t131 * ((ddreal(56) * s25 + t149) * s25 - ddreal(9) 
     &* t9) + ((s16 * (((ddreal(97) * s16 + t161) * s25 + t234) * s25 + 
     &t13 * t165) + (t177 + ((ddreal(80) * s16 + t163) * s25 + t32) * s2
     &5 + t170) * s26) * s26 + t9 * (((t98 + t190) * s25 + ddreal(98) * 
     &t9) * s25 - t198)) * s26 - t11) * s34
      t84 = t47 * (-t42 * t81 * t4 + (((((ddreal(77) * s16 + t171) * s25
     & + t192) * s26 + ((ddreal(216) * s16 + t188) * s25 + ddreal(166) *
     & t9) * s25 + t244) * s26 + s16 * (((ddreal(212) * s16 + t179) * s2
     &5 + ddreal(97) * t9) * s25 - t185)) * s26 - t9 * (-t110 + s16) * (
     &(ddreal(64) * s16 + t84) * s25 + t46)) * s26 + t197 * t64 - t138 *
     & t1) * t152
      t7 = s56 * (s56 * (s34 * (s34 * (-t197 * s34 + t172) - t148 * t8) 
     &+ t1 * t194 + t26 * t206) + t6 * (s16 * (t240 + t10) + t153 * (t5 
     &* t7 - s34 + t143)) - (ddreal(24) * s26 * t148 + s34 * (s34 * (t25
     &6 * t270 + t207 + ddreal(60) * t247) - s25 * (s25 * (t235 + t70) +
     & ddreal(78) * t9) - t170 - t105 * t172) + t232 + t33) * s34 + t26 
     &* (s25 * (s25 * (s16 * t29 - t163) - t76) - t96)) + s16 * (s25 * (
     &s25 * (s25 * (-s25 * t93 + t23 * t9) - ddreal(57) * t13) + t162) +
     & t262 * t178 * t111 + t257 * t49) + t1 * (t156 + t11) + (s34 * (s3
     &4 * (s34 * (ddreal(45) * t62 + ddreal(150) * t151 + t166 + t112) -
     & s26 * (t241 + t228) - t244 - t89) + s26 * (t119 * t24 + t130) + t
     &133 + t266) + s26 * (-t238 * t7 - t258) - t12 + t26 * (s25 * (-t18
     &7 + t225) - t175)) * s34 - t251 * (t14 + t128 + t142)
      t7 = s56 * t7 - t139 * t67 * t103 + s16 * (t239 * s26 + t82 * (t17
     & * (-ddreal(68) * s16 + ddreal(92) * s25) + t94)) - ddreal(20) * s
     &16 * (-t257 * (s26 + t111) + t20) + ((((s34 * (s34 * (t141 - t54 +
     & t14) + s26 * (t111 * t158 - t252) + t189 + t214) + ((t253 + t252)
     & * s26 + t248) * s26 + t203 - t254) * s34 + s16 * (((s25 * t224 - 
     &ddreal(378) * s16) * s25 - ddreal(121) * t9) * s25 + t185) + (-t17
     &8 * t227 - t230) * s26 - ddreal(40) * t193 * t65 - ddreal(40) * t1
     &8) * s34 + s16 * ((((ddreal(404) * s16 + ddreal(314) * s25) * s25 
     &+ ddreal(325) * t9) * s25 + ddreal(140) * t13) * s25 - t10 * t195)
     & + ddreal(20) * t62 * (t24 + t177) + ddreal(20) * t21 + t8 * t267 
     &* s26) * s34 + t62 * (-ddreal(80) * s26 * t148 - t119 * t238) - t8
     & * ((t12 + t26 * ((t187 - t225) * s25 + t175)) * s26 - t11 + t20) 
     &+ t26 * ((((-t221 - t78) * s25 + ddreal(236) * t9) * s25 - t180) *
     & s25 - ddreal(67) * t10)) * s34
      t7 = s56 * t7 + s16 * (t17 * (((((-t261 - t122) * s25 + ddreal(63)
     & * t9) * s25 - t45) * s25 - t237) * s25 + t255) + t205 * t115 * t1
     &47 + ddreal(20) * t262 * t200 + t262 * t49 * t63 + t69 * t141 * t1
     &15 * t116) + (t100 * s34 - t47 * (s16 * (((((t123 + t191) * s25 - 
     &ddreal(144) * t9) * s25 + t13 * t158) * s25 + ddreal(129) * t10) *
     & s25 - t68) + ((t119 * ((((t98 + t43) * s25 - t108) * s25 - t40) *
     & s25 + t33) + t258) * s26 + t8 * (((((t73 + t2) * s25 + t229) * s2
     &5 - ddreal(141) * t13) * s25 + ddreal(93) * t10) * s25 + t255)) * 
     &s26)) * s34
      t6 = s56 * t7 + s16 * (t134 * t93 + t239 * t65) - t119 * s16 * (s1
     &6 * t18 * t19 - t262 * t263) + t31 * t256 * (t115 * t16 * t164 + t
     &155) + t50 * t26 * (t116 * t257 + t184) + t6 * s16 * (t20 * t9 + t
     &262 * t64) + s34 * (t72 * s34 - t199 * (-t11 * t50 + t62 * ((((t75
     & + t159) * s25 - ddreal(86) * t9) * s25 - ddreal(99) * t13) * s25 
     &+ ddreal(150) * t10) + t211 + t173 * (t17 * (ddreal(19) * t17 - dd
     &real(155) * t9) + ddreal(74) * t136 + t183) + ddreal(24) * t150 + 
     &ddreal(105) * t138 + t77 * t29 - ddreal(122) * t209))
      t7 = t26 * (t212 * t116 * t65 + t236 * (s25 * t206 + t25))
      t12 = ddreal(1) / s16
      t14 = ddreal(1) / (t35)**2
      t20 = ddreal(1) / (t39)**2
      t23 = ddreal(1) / (t34)**2
      t3 = t61 * (t99 + (-t157 * (t31 * (t82 + t231) + t233) + t60 * (((
     &t215 + t59) * s25 - t44) * s25 + t40)) * s26 + t38 * s34 + s56 * t
     &97 + t160 * t28 + t167 * t57 + t169 * ((((((t55 + t51) * s25 - ddr
     &eal(81) * t9) * s25 + ddreal(98) * t13) * s25 - ddreal(43) * t10) 
     &* s25 - t269 * t3 * t119) * s25 + t1 * t11 + t31 * ((-s26 * t31 * 
     &t71 + ddreal(20) * t132 * t17 * t62 - s25 * t124 - t153 * t226 - t
     &18 * t51 - ddreal(19) * t21 + ddreal(51) * t236 - t268 + ddreal(17
     &) * t81) * s34 + t138) - t8 * (t154 * (t103 * t50 + t261 * t62 + d
     &dreal(14) * t144 + t213 + ddreal(20) * t82) - t182 * t92) + t152 *
     & (t259 * s26 + t152 * (s16 * (t249 + t122) + (s26 * t42 + t210 + t
     &219) * s25) + t176 * t62 + t186 * t48 + ddreal(59) * t131 + ddreal
     &(62) * t18 + t237 + t87)) + t62 * (t69 * t85 + t79) + t95 + t264 *
     & t92 * t63 + t264 * t91 * t65 - t69 * t3 * t64) * MBL1011(1)
      t2 = (t31 * t240 * s16 * (t186 * t101 + t263) + s34 * (-t47 * t199
     & * (s16 * (t62 * (((t123 - t120) * s25 - t140) * s25 + t135) + t23
     &6 * t58 + t81 * t158 + t218 * t10 + t173 * (((t202 - t114) * s25 +
     & t222) * s25 + t83)) + t8 * (t150 - t11)) + t154 * (((t65 * (-t111
     & * t50 - t105) + t82 * (t2 + s16) * t139) * s26 + t62 * (s26 * t21
     &6 + t217) + t81) * s34 - t1 * (-t62 * t63 + t209) + (((s26 * t208 
     &+ t204) * s26 + s16 * (((t98 + t242) * s25 - t201) * s25 - t96)) *
     & s26 - t82 * (t50 * t66 + ddreal(65) * t26)) * s26 - t211 + t30 * 
     &t5 * t64) + t32 - t84 + s26 * (t1 * t65 - t86 + s26 * (s26 * t220 
     &- t223)) * t157 + t79 * t62) + s56 * t6 + t257 * s16 * (t16 * t17 
     &+ t62 * (t62 + t147)) + t7) * MBL1110(1) + t27 + t3
      t2 = (t1 * s25 * (-t168 * (t1 * t17 - t8 * s25 * (s16 + s34) + t15
     &2 + t9 - s34 * t128) + t184 + t4 * t62 * (-t245 + t181)) - ddreal(
     &13) * t77 * t125 + t22 + t198 * t21) * MBL1001(0) + t2 * t12
      t3 = (t61)**2
      t3 = t12 * (t20 * (-t61 * t3 * t113 * MBL1111D2(0) + (t21 * t52 - 
     &t107) * MBL1101(0) * t23) + (t113)**2 * MBL1110(0) * t14)
      result = t23 * t20 * ((t1 * ((t82 * (-t37 + t110) - t8 * t5 * t154
     & - ((-t152 - t126) * s34 + t31 * t104) * s34 - t19 * t4 + t10) * t
     &169 + t174 + t4 * t53 * t65) + ((((s26 * ((t145 + t15 + t143) * s2
     &6 + t146) - s34 * t62) * s34 - t62 * (t137 * t31 + t141 * t5) - t6
     &7 - t121 * t129 * s26) * s34 + t1 * (t126 * t65 + t81) + s26 * (s2
     &6 * t118 + t127) + t46 * t19) * s34 - t62 * (t104 * t105 + t102) -
     & t106 * t17 - ddreal(7) * t77 - t36 * t26 * ((s25 * t88 - t9) * s2
     &5 + t90)) * s34 + t109 * s56 + t62 * (t80 + t74 + t56) + t184 + t1
     &17 * (((-t196 - t43) * s25 - t46) * s25 + t52) * s26 + t134 - t41 
     &* t19) * I300s16(1) + t14 * t2) / ddreal(12) + t3 / ddreal(6)

           intHLs16s25s26s34s56x1411D6eps0 = result
       end function intHLs16s25s26s34s56x1411D6eps0

       function intHLs16s25s26s34s56x1411D6eps1()
           implicit none
           type(dd_complex) :: intHLs16s25s26s34s56x1411D6eps1
           type(dd_complex) ::  t1,t2,t3,t4

           type(dd_complex) :: result

      t1 = s25 + s26 + s56
      t2 = ddreal(3)
      t3 = s25 + s26 - s34 + s56
      t4 = ddreal(1) / ddreal(12)
      t3 = ddreal(1) / (t3)**2
      result = ddreal(1) / s16 * t3 * (t4 * (t1 * t2 - s34) + (t1)**2 * 
     &MBL1110(1) / ddreal(6))

           intHLs16s25s26s34s56x1411D6eps1 = result
       end function intHLs16s25s26s34s56x1411D6eps1

       function intHLs250000x0111D2eps0()
           implicit none
           type(dd_complex) :: intHLs250000x0111D2eps0

           type(dd_complex) :: result

      result = ddreal(3) / ddreal(2) + s25 * I300s25(1) / ddreal(2)

           intHLs250000x0111D2eps0 = result
       end function intHLs250000x0111D2eps0

       function intHLs250000x0111D2eps1()
           implicit none
           type(dd_complex) :: intHLs250000x0111D2eps1

           type(dd_complex) :: result

      result = ddreal(1) / ddreal(2)

           intHLs250000x0111D2eps1 = result
       end function intHLs250000x0111D2eps1

       function intHLs250000x0112D2eps0()
           implicit none
           type(dd_complex) :: intHLs250000x0112D2eps0

           type(dd_complex) :: result

      result = -ddreal(2) / s25 - I300s25(1)

           intHLs250000x0112D2eps0 = result
       end function intHLs250000x0112D2eps0

       function intHLs250000x0112D2eps1()
           implicit none
           type(dd_complex) :: intHLs250000x0112D2eps1

           type(dd_complex) :: result

      result = -ddreal(1) / s25

           intHLs250000x0112D2eps1 = result
       end function intHLs250000x0112D2eps1

       function intHLs250000x0113D4eps0()
           implicit none
           type(dd_complex) :: intHLs250000x0113D4eps0

           type(dd_complex) :: result

      result = -ddreal(1) / s25 / ddreal(2) - I300s25(1) / ddreal(4)

           intHLs250000x0113D4eps0 = result
       end function intHLs250000x0113D4eps0

       function intHLs250000x0113D4eps1()
           implicit none
           type(dd_complex) :: intHLs250000x0113D4eps1

           type(dd_complex) :: result

      result = -ddreal(1) / s25 / ddreal(4)

           intHLs250000x0113D4eps1 = result
       end function intHLs250000x0113D4eps1

       function intHLs250000x0121D2eps0()
           implicit none
           type(dd_complex) :: intHLs250000x0121D2eps0

           type(dd_complex) :: result

      result = ddreal(4) / s25 + I300s25(0) + ddreal(2) * I300s25(1)

           intHLs250000x0121D2eps0 = result
       end function intHLs250000x0121D2eps0

       function intHLs250000x0121D2eps1()
           implicit none
           type(dd_complex) :: intHLs250000x0121D2eps1

           type(dd_complex) :: result

      result = ddreal(2) / s25 + I300s25(1)

           intHLs250000x0121D2eps1 = result
       end function intHLs250000x0121D2eps1

       function intHLs250000x0122D4eps0()
           implicit none
           type(dd_complex) :: intHLs250000x0122D4eps0

           type(dd_complex) :: result

      result = -ddreal(3) / ddreal(2) / s25 - I300s25(1) / ddreal(2)

           intHLs250000x0122D4eps0 = result
       end function intHLs250000x0122D4eps0

       function intHLs250000x0122D4eps1()
           implicit none
           type(dd_complex) :: intHLs250000x0122D4eps1

           type(dd_complex) :: result

      result = -ddreal(1) / s25 / ddreal(2)

           intHLs250000x0122D4eps1 = result
       end function intHLs250000x0122D4eps1

       function intHs160000x0111D0eps0()
           implicit none
           type(dd_complex) :: intHs160000x0111D0eps0

           type(dd_complex) :: result

      result = I300s16(0)

           intHs160000x0111D0eps0 = result
       end function intHs160000x0111D0eps0

       function intHs160000x0111D0eps1()
           implicit none
           type(dd_complex) :: intHs160000x0111D0eps1

           type(dd_complex) :: result

      result = I300s16(1)

           intHs160000x0111D0eps1 = result
       end function intHs160000x0111D0eps1

       function intHs160000x0112D2eps0()
           implicit none
           type(dd_complex) :: intHs160000x0112D2eps0

           type(dd_complex) :: result

      result = -ddreal(2) / s16 - I300s16(1)

           intHs160000x0112D2eps0 = result
       end function intHs160000x0112D2eps0

       function intHs160000x0112D2eps1()
           implicit none
           type(dd_complex) :: intHs160000x0112D2eps1

           type(dd_complex) :: result

      result = -ddreal(1) / s16

           intHs160000x0112D2eps1 = result
       end function intHs160000x0112D2eps1

       function intHs160000x0121D2eps0()
           implicit none
           type(dd_complex) :: intHs160000x0121D2eps0

           type(dd_complex) :: result

      result = ddreal(4) / s16 + I300s16(0) + ddreal(2) * I300s16(1)

           intHs160000x0121D2eps0 = result
       end function intHs160000x0121D2eps0

       function intHs160000x0121D2eps1()
           implicit none
           type(dd_complex) :: intHs160000x0121D2eps1

           type(dd_complex) :: result

      result = ddreal(2) / s16 + I300s16(1)

           intHs160000x0121D2eps1 = result
       end function intHs160000x0121D2eps1

       function intHs160000x0211D2eps0()
           implicit none
           type(dd_complex) :: intHs160000x0211D2eps0

           type(dd_complex) :: result

      result = -ddreal(2) / s16 - I300s16(1)

           intHs160000x0211D2eps0 = result
       end function intHs160000x0211D2eps0

       function intHs160000x0211D2eps1()
           implicit none
           type(dd_complex) :: intHs160000x0211D2eps1

           type(dd_complex) :: result

      result = -ddreal(1) / s16

           intHs160000x0211D2eps1 = result
       end function intHs160000x0211D2eps1

       function intHs160s26s34s56x1011D2eps0()
           implicit none
           type(dd_complex) :: intHs160s26s34s56x1011D2eps0
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = ddreal(1) / ddreal(2)
      result = t1 * (-(s16 + s26 - s34 + s56) * MB1011(1) + MB1001(0) + 
     &ddreal(1))

           intHs160s26s34s56x1011D2eps0 = result
       end function intHs160s26s34s56x1011D2eps0

       function intHs160s26s34s56x1011D2eps1()
           implicit none
           type(dd_complex) :: intHs160s26s34s56x1011D2eps1

           type(dd_complex) :: result

      result = ddreal(1) / ddreal(2)

           intHs160s26s34s56x1011D2eps1 = result
       end function intHs160s26s34s56x1011D2eps1

       function intHs160s26s34s56x1012D2eps0()
           implicit none
           type(dd_complex) :: intHs160s26s34s56x1012D2eps0
           type(dd_complex) ::  t1,t2

           type(dd_complex) :: result

      t1 = s16 + s26 + s56
      t2 = s16 + s26 - s34 + s56
      t1 = ddreal(1) / t1
      result = ddreal(2) * t2 * MB1011(1) * t1 + t1 * (t2 * MB1011(0) - 
     &MB1001(0))

           intHs160s26s34s56x1012D2eps0 = result
       end function intHs160s26s34s56x1012D2eps0

       function intHs160s26s34s56x1012D2eps1()
           implicit none
           type(dd_complex) :: intHs160s26s34s56x1012D2eps1
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s16 + s26 + s56
      t1 = ddreal(1) / t1
      result = t1 * ((s16 + s26 - s34 + s56) * MB1011(1) - ddreal(1))

           intHs160s26s34s56x1012D2eps1 = result
       end function intHs160s26s34s56x1012D2eps1

       function intHs160s26s34s56x1013D4eps0()
           implicit none
           type(dd_complex) :: intHs160s26s34s56x1013D4eps0
           type(dd_complex) ::  t1,t2,t3

           type(dd_complex) :: result

      t1 = s16 + s26 - s34 + s56
      t2 = s16 + s26 + s56
      t3 = ddreal(1) / (t2)**2
      result = t3 * (-(-ddreal(2) * s34 + ddreal(3) * t2) * MB1001(0) / 
     &ddreal(4) + ddreal(3) / ddreal(2) * (t1)**2 * MB1011(1)) - t1 * t3
     & * (-t1 * MB1011(0) + ddreal(1)) / ddreal(2)

           intHs160s26s34s56x1013D4eps0 = result
       end function intHs160s26s34s56x1013D4eps0

       function intHs160s26s34s56x1013D4eps1()
           implicit none
           type(dd_complex) :: intHs160s26s34s56x1013D4eps1
           type(dd_complex) ::  t1,t2,t3

           type(dd_complex) :: result

      t1 = s16 + s26 + s56
      t2 = s16 + s26 - s34 + s56
      t3 = ddreal(1) / (t1)**2
      result = t3 * (s34 / ddreal(2) - ddreal(3) / ddreal(4) * t1 + (t2)
     &**2 * MB1011(1) / ddreal(2))

           intHs160s26s34s56x1013D4eps1 = result
       end function intHs160s26s34s56x1013D4eps1

       function intHs160s26s34s56x1020D2eps0()
           implicit none
           type(dd_complex) :: intHs160s26s34s56x1020D2eps0
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = ddreal(1) / ddreal(2)
      result = t1 * (-(s16 + s26 + s56) * MB1011(1) + MB1001(0))

           intHs160s26s34s56x1020D2eps0 = result
       end function intHs160s26s34s56x1020D2eps0

       function intHs160s26s34s56x1020D2eps1()
           implicit none
           type(dd_complex) :: intHs160s26s34s56x1020D2eps1

           type(dd_complex) :: result

      result = ddreal(1) / ddreal(2)

           intHs160s26s34s56x1020D2eps1 = result
       end function intHs160s26s34s56x1020D2eps1

       function intHs160s26s34s56x1021D2eps0()
           implicit none
           type(dd_complex) :: intHs160s26s34s56x1021D2eps0
           type(dd_complex) ::  t1,t2

           type(dd_complex) :: result

      t1 = s16 + s26 + s56
      t2 = ddreal(1) / t1
      result = (s34 * MB1011(0) - (-ddreal(2) * s34 + t1) * MB1011(1) + 
     &MB1001(0)) * t2

           intHs160s26s34s56x1021D2eps0 = result
       end function intHs160s26s34s56x1021D2eps0

       function intHs160s26s34s56x1021D2eps1()
           implicit none
           type(dd_complex) :: intHs160s26s34s56x1021D2eps1
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s16 + s26 + s56
      t1 = ddreal(1) / t1
      result = t1 * (s34 * MB1011(1) + ddreal(1))

           intHs160s26s34s56x1021D2eps1 = result
       end function intHs160s26s34s56x1021D2eps1

       function intHs160s26s34s56x1022D4eps0()
           implicit none
           type(dd_complex) :: intHs160s26s34s56x1022D4eps0
           type(dd_complex) ::  t1,t2,t3

           type(dd_complex) :: result

      t1 = s16 + s26 + s56
      t2 = s16 + s26 - s34 + s56
      t3 = ddreal(1) / (t1)**2
      result = -t3 * (-(MB1001(0) + ddreal(1)) * (-ddreal(2) * s34 + t1)
     & + (-ddreal(6) * s34 + t1) * t2 * MB1011(1)) / ddreal(2) + s34 * t
     &2 * MB1011(0) * t3

           intHs160s26s34s56x1022D4eps0 = result
       end function intHs160s26s34s56x1022D4eps0

       function intHs160s26s34s56x1022D4eps1()
           implicit none
           type(dd_complex) :: intHs160s26s34s56x1022D4eps1
           type(dd_complex) ::  t1,t2

           type(dd_complex) :: result

      t1 = s16 + s26 + s56
      t2 = ddreal(1) / (t1)**2
      result = t2 * (-s34 + t1 / ddreal(2) + s34 * (s16 + s26 - s34 + s5
     &6) * MB1011(1))

           intHs160s26s34s56x1022D4eps1 = result
       end function intHs160s26s34s56x1022D4eps1

       function intHs160s26s34s56x1031D4eps0()
           implicit none
           type(dd_complex) :: intHs160s26s34s56x1031D4eps0
           type(dd_complex) ::  t1,t2

           type(dd_complex) :: result

      t1 = s16 + s26 + s56
      t2 = ddreal(1) / (t1)**2
      result = s34 * t2 * (s34 * MB1011(0) + ddreal(1)) / ddreal(2) + t2
     & * ((ddreal(2) * s34 + t1) * MB1001(0) - ((s16)**2 + (s26)**2 + (s
     &56)**2 - ddreal(6) * (s34)**2 + ddreal(2) * (s16 + s26 + s34) * s5
     &6 + ddreal(2) * (s16 + s26) * s34 + ddreal(2) * s16 * s26) * MB101
     &1(1)) / ddreal(4)

           intHs160s26s34s56x1031D4eps0 = result
       end function intHs160s26s34s56x1031D4eps0

       function intHs160s26s34s56x1031D4eps1()
           implicit none
           type(dd_complex) :: intHs160s26s34s56x1031D4eps1
           type(dd_complex) ::  t1,t2

           type(dd_complex) :: result

      t1 = s16 + s26 + s56
      t2 = ddreal(1) / (t1)**2
      result = t2 * (s34 / ddreal(2) + t1 / ddreal(4) + (s34)**2 * MB101
     &1(1) / ddreal(2))

           intHs160s26s34s56x1031D4eps1 = result
       end function intHs160s26s34s56x1031D4eps1

       function intHs16s25s26s34s56x1110D2eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1110D2eps0
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = ddreal(1) / ddreal(2)
      result = t1 * (s25 * MB1110(1) - (s16 + s26 + s56) * MB1011(1) + M
     &B1001(0) + ddreal(1))

           intHs16s25s26s34s56x1110D2eps0 = result
       end function intHs16s25s26s34s56x1110D2eps0

       function intHs16s25s26s34s56x1110D2eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1110D2eps1

           type(dd_complex) :: result

      result = ddreal(1) / ddreal(2)

           intHs16s25s26s34s56x1110D2eps1 = result
       end function intHs16s25s26s34s56x1110D2eps1

       function intHs16s25s26s34s56x1111D2eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1111D2eps0
           type(dd_complex) ::  t1,t2,t3,t4,t5,t6,t7,t8

           type(dd_complex) :: result

      t1 = s16 + s26 - s34 + s56
      t2 = s12 + s15 + s25
      t3 = ddreal(2)
      t4 = s16 * s25 + (-s34 + s16 + s25 + s26) * s26 + (s26 * t3 + s16 
     &+ s25 - s34 + s56) * s56
      t5 = s16 + s26 + s56
      t6 = s16 - s25
      t7 = s16 + s25 + s26 - s34 + s56
      t8 = ddreal(1) / t4
      t4 = ddreal(1) / t4
      t4 = t4 * t1
      t2 = ddreal(1) / t2
      result = -t8 * (t1 * ((t1 * MB1111(0) + I300s16(0)) * s16 - t5 * M
     &B1011(0) - t7 * MB1110(0)) + (s16 * t6 + s26 * t6 + (-s16 * t3 - s
     &25 - s26 + s34) * s34 + (s16 - s34 - s25) * s56) * MB1101(0)) / dd
     &real(2) - t4 * ((t1 * MB1111(1) + I300s16(1)) * s16 - t5 * MB1011(
     &1) - t7 * MB1110(1)) - t4 * t3 * (s12 + s15 + s16 + s25 + s26 - s3
     &4 + s56) * t2

           intHs16s25s26s34s56x1111D2eps0 = result
       end function intHs16s25s26s34s56x1111D2eps0

       function intHs16s25s26s34s56x1111D2eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1111D2eps1
           type(dd_complex) ::  t1,t2,t3,t4

           type(dd_complex) :: result

      t1 = s16 + s26 - s34 + s56
      t2 = s12 + s15 + s25
      t3 = ddreal(2)
      t3 = s16 * s25 + (-s34 + s16 + s25 + s26) * s26 + (s26 * t3 + s16 
     &+ s25 - s34 + s56) * s56
      t4 = ddreal(1) / t3
      t3 = ddreal(1) / t3
      t2 = ddreal(1) / t2
      result = -t4 * t1 * ((t1 * MB1111(1) + I300s16(1)) * s16 - (s16 + 
     &s26 + s56) * MB1011(1) - (s16 + s25 + s26 - s34 + s56) * MB1110(1)
     &) / ddreal(2) - t1 * (s12 + s15 + s16 + s25 + s26 - s34 + s56) * t
     &2 * t3

           intHs16s25s26s34s56x1111D2eps1 = result
       end function intHs16s25s26s34s56x1111D2eps1

       function intHs16s25s26s34s56x1112D2eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1112D2eps0
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s16 + s26 - s34 + s56
      t1 = ddreal(1) / t1
      result = ddreal(1) / s16 * ((s25 * MB1101(0) - (s16 + s25 + s26 - 
     &s34 + s56) * MB1111D2(0)) * t1 + MB1011(0))

           intHs16s25s26s34s56x1112D2eps0 = result
       end function intHs16s25s26s34s56x1112D2eps0

       function intHs16s25s26s34s56x1112D2eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1112D2eps1

           type(dd_complex) :: result

      result = MB1011(1) / s16

           intHs16s25s26s34s56x1112D2eps1 = result
       end function intHs16s25s26s34s56x1112D2eps1

       function intHs16s25s26s34s56x1112D4eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1112D4eps0
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           type(dd_complex) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t4
           type(dd_complex) ::  t5,t6,t7,t8,t9

           type(dd_complex) :: result

      t1 = s16 - s25
      t2 = ddreal(2)
      t3 = t2 * s16
      t4 = s25 + t3
      t5 = s16 + s25
      t6 = s25 * t5
      t7 = (s16)**2
      t8 = ddreal(3) * s16
      t9 = (s34)**2
      t10 = s34 * t9
      t11 = s16 * s25
      t12 = t1 * s16
      t13 = t1 * s26
      t14 = (-t13 - t12) * s16
      t15 = (s26 * t4 + (-s25 * t2 - s26 + s34 - t8) * s34 + t6 + ddreal
     &(3) * t7) * s34
      t16 = -(-s25 * s34 - t3 * s34 - t11 + t7 + t9) * s56 + t14 + t15
      t17 = (s25)**2
      t18 = -t2 * (s34 * t5 + t11) + t17 + t7 + t9
      t19 = s16 + s25 - s34
      t20 = t2 * s26
      t11 = (-s34 + t5 + s26) * s26 + (t20 + t19 + s56) * s56 + t11
      t21 = s16 + s26
      t22 = s25 + s26
      t23 = ddreal(4) * s16
      t24 = (s26)**2
      t25 = t3 + t22
      t26 = t1 * t7
      t27 = s16 + s25 + s26 - s34 + s56
      t18 = ddreal(1) / t18
      t28 = ddreal(1) / t11
      t11 = ddreal(1) / t11
      t29 = ddreal(1) / ddreal(2)
      result = t29 * ((s16 + s26 - s34 + s56) * t27 * MB1111D2(0) * t28 
     &+ ((s34 * (s26 * t5 + (-t2 * t5 - s26 + s34) * s34 + s56 * t19 + t
     &17 + t7) * MB1101(0) + t27 * ((t25 - s34) * s34 - (s16 - s25 - s34
     &) * s56 - t12 - t13) * MB1110(1)) * s25 + (((t4 - s34) * s34 - t12
     &) * s56 + t14 + t15) * MB1001(0) + (t2 * (t26 * s26 - t10 * t25) -
     & (-t1 * t24 - t26) * s16 - ((-t9 - ddreal(6) * s16 * t21 - (ddreal
     &(5) * s26 + t8) * s25 - t17 - t24) * s34 + t21 * ((t23 + t20) * s1
     &6 + ddreal(3) * s25 * t22)) * s34 - (-(-(ddreal(3) * s25 + t3) * s
     &34 + t12 + t9) * s56 + t2 * (t14 + t10) + s34 * (ddreal(6) * s25 *
     & s26 + t23 * s26 - (ddreal(6) * s16 + ddreal(5) * s25 + t20) * s34
     & + ddreal(3) * t6 + ddreal(6) * t7)) * s56) * MB1011(1) - s16 * t1
     &6 * I300s16(1)) * t11 * t18) - t16 * t18 * t11

           intHs16s25s26s34s56x1112D4eps0 = result
       end function intHs16s25s26s34s56x1112D4eps0

       function intHs16s25s26s34s56x1112D4eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1112D4eps1

           type(dd_complex) :: result

      result = ddreal(0)

           intHs16s25s26s34s56x1112D4eps1 = result
       end function intHs16s25s26s34s56x1112D4eps1

       function intHs16s25s26s34s56x1113D4eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1113D4eps0
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           type(dd_complex) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           type(dd_complex) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           type(dd_complex) ::  t42,t43,t5,t6,t7,t8,t9

           type(dd_complex) :: result

      t1 = (s25)**2
      t2 = s25 * t1
      t3 = (s16)**2
      t4 = s16 * t3
      t5 = t1 + t3
      t6 = s16 + s25
      t7 = s25 * t6
      t8 = ddreal(3)
      t9 = ddreal(2)
      t10 = s34 * t6
      t11 = (s34)**2
      t12 = t1 + t3 + t11
      t13 = s26 * t5 + ((t6 * t8 + s26 - s34) * s34 - t8 * (t7 + t3) - t
     &9 * s26 * t6) * s34 + (-t10 * t9 + t12) * s56 + t2 + t4
      t14 = s16 * s25
      t10 = -t9 * (t10 + t14) + t12
      t12 = s16 + s25 - s34
      t15 = (-s34 + t6 + s26) * s26 + (s26 * t9 + s56 + t12) * s56 + t14
      t16 = s16 * t6
      t17 = t9 * t1
      t18 = t9 * s16
      t19 = t9 * s25
      t20 = ddreal(4)
      t21 = t20 * s26
      t22 = (s26)**2
      t23 = -s16 + s25
      t24 = s16 + s26 + s56
      t25 = s16 * t23
      t26 = ddreal(6)
      t27 = t8 * s25
      t28 = (-s16 * t20 - t27) * s25 + t3
      t29 = s16 + s26
      t30 = (s25 + t18) * s26
      t31 = ddreal(5) * s16
      t32 = t20 * s25
      t33 = t8 * t1
      t34 = t8 * s26
      t35 = t20 * t3
      t36 = (s16 - t27) * s26
      t37 = s16 * t26
      t37 = (-(-(s34 * t9 + t27) * s16 - (s25 - s34) * s34 + t3) * s56 +
     & t9 * s34 * t11 + t8 * s16 * (-t36 + t1) - ((t34 + t31 + t32) * s3
     &4 - (t37 + t27) * s26 + t20 * t25 - t17) * s34 - t4 + t35 * s25) *
     & s56 - s16 * ((t28 * t9 + t36 * t8) * s26 - t37 * t1) - (((-t27 - 
     &t21 - t18 + s34) * s34 + (ddreal(10) * s16 + ddreal(8) * s25) * s2
     &6 + t8 * (t1 + t22)) * s34 + (-t20 * (-t25 * t9 + t1) - t30 * t8) 
     &* s26 + t9 * s16 * ((s16 + s25) * s25 + t3) - t2) * s34 - t23 * t4
      t38 = s16 + t19
      t39 = t23 * s25
      t40 = s16 + s25 + s26 - s34 + s56
      t24 = ddreal(1) / t24
      t10 = ddreal(1) / t10
      t15 = ddreal(1) / t15
      t41 = ddreal(1) / s16
      t42 = t10 * t15
      t10 = t13 * I300s16(1) * t10 + (t40)**2 * MB1111D2(0) * t41
      t43 = t42 * t41
      t6 = t43 * (((-s26 * t23 - t39) * s25 + (s26 * t38 + (-s26 - t27 -
     & t18 + s34) * s34 + t16 + t33) * s34 + ((t38 - s34) * s34 - t39) *
     & s56) * MB1101(0) - t40 * (s26 * t6 + (-t6 * t9 - s26 + s34) * s34
     & + s56 * t12 + t1 + t3) * MB1110(1)) * s25
      t12 = -ddreal(1) / ddreal(2)
      result = t12 * (t10 * t15 + t24 * (t42 * (-(t9 * ((s56)**2 + t14 +
     & t22) * s25 + (t17 + t16) * s26 - ((-s16 * t8 - s26 + s34 - t19) *
     & s34 + t8 * (s25 * s26 + t3) + t7 + t18 * s26) * s34 + ((-s34 * t8
     & + t21) * s25 - t9 * (s16 * s34 - t1) + t11 + t16) * s56 - t3 * t2
     &3) * MB1001(0) - (t4 * s25 * (-t27 + t18) - s16 * ((s16 * (t1 * t2
     &6 - t25) + (-t36 - t28) * s26) * s26 + (t3)**2) - ((((s16 - s26) *
     & s34 + t9 * s26 * t29 - (s16 - t34) * s25 - t35) * s34 + t26 * t4 
     &- s16 * t1 - t22 * (t31 + t32 + s26) - t33 * s26) * s34 - t29 * (t
     &20 * s16 * t5 - ((s16 - t19) * (-s25 + t18) + t30) * s26 - t2 - t2
     &7 * t3)) * s34 - t37 * s56) * MB1011(1) * t41) - (s16 + s26 - s34 
     &+ s56) * MB1011(0) * t41) + t6) - t43 * t13

           intHs16s25s26s34s56x1113D4eps0 = result
       end function intHs16s25s26s34s56x1113D4eps0

       function intHs16s25s26s34s56x1113D4eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1113D4eps1
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s16 + s26 + s56
      t1 = ddreal(1) / t1
      result = ddreal(1) / s16 * t1 * ((s16 + s26 - s34 + s56) * MB1011(
     &1) - ddreal(1)) / ddreal(2)

           intHs16s25s26s34s56x1113D4eps1 = result
       end function intHs16s25s26s34s56x1113D4eps1

       function intHs16s25s26s34s56x1114D6eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1114D6eps0
           type(dd_complex) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           type(dd_complex) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           type(dd_complex) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           type(dd_complex) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           type(dd_complex) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           type(dd_complex) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           type(dd_complex) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           type(dd_complex) ::  t175,t176,t177,t178,t179,t18,t180,t181,t182,t183,t184,t185
           type(dd_complex) ::  t186,t187,t188,t189,t19,t190,t191,t192,t193,t194,t195,t196
           type(dd_complex) ::  t197,t198,t199,t2,t20,t200,t201,t202,t203,t204,t205,t206
           type(dd_complex) ::  t207,t208,t209,t21,t210,t211,t212,t213,t214,t215,t216,t217
           type(dd_complex) ::  t218,t219,t22,t220,t221,t222,t223,t224,t225,t226,t227,t228
           type(dd_complex) ::  t229,t23,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239
           type(dd_complex) ::  t24,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249,t25
           type(dd_complex) ::  t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t26,t260
           type(dd_complex) ::  t261,t262,t263,t264,t265,t266,t267,t268,t269,t27,t270,t271
           type(dd_complex) ::  t272,t273,t274,t275,t276,t277,t278,t279,t28,t280,t281,t282
           type(dd_complex) ::  t283,t284,t285,t286,t287,t288,t289,t29,t290,t291,t292,t293
           type(dd_complex) ::  t294,t295,t296,t3,t30,t31,t32,t33,t34,t35,t36,t37
           type(dd_complex) ::  t38,t39,t4,t40,t41,t42,t43,t44,t45,t46,t47,t48
           type(dd_complex) ::  t49,t5,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59
           type(dd_complex) ::  t6,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t7
           type(dd_complex) ::  t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t8,t80
           type(dd_complex) ::  t81,t82,t83,t84,t85,t86,t87,t88,t89,t9,t90,t91
           type(dd_complex) ::  t92,t93,t94,t95,t96,t97,t98,t99

           type(dd_complex) :: result

      t1 = s16 - s25
      t2 = s16 + s25
      t3 = (s16)**2
      t4 = (t3)**2
      t5 = (t4)**2
      t6 = s16 * t5
      t7 = s16 * t3
      t8 = t7 * t4
      t9 = t3 * t4
      t10 = s16 * t4
      t11 = t3 * s25
      t12 = t11 * t2
      t13 = ddreal(8)
      t14 = ddreal(6)
      t15 = ddreal(5)
      t16 = (s25)**2
      t17 = (t16)**2
      t18 = s25 * t16
      t19 = t16 * t17
      t20 = s25 * t17
      t21 = t15 * s16
      t22 = s25 * t14
      t23 = ddreal(7) * t4
      t24 = ddreal(3)
      t25 = ddreal(18) * s25
      t26 = ddreal(27) * s16
      t27 = t24 * t3
      t28 = ddreal(46) * t10
      t29 = ((((-t25 - t26) * s25 - t27) * s25 - ddreal(87) * t7) * s25 
     &- ddreal(31) * t4) * s25 + t28
      t30 = ddreal(52)
      t31 = ddreal(33)
      t32 = t3 * t15
      t33 = ddreal(65) * t4
      t34 = ddreal(2)
      t35 = ddreal(7) * t7
      t36 = ((t18 * t34 + t35) * s25 - ddreal(14) * t4) * s25 + t10 * t1
     &5
      t37 = t24 * s16
      t38 = t37 - t22
      t39 = ddreal(16) * t3
      t40 = (s25 * t38 - t39) * s25 + t35
      t41 = s16 + s26
      t42 = ddreal(47) * s25
      t43 = ddreal(24) * s16
      t44 = ddreal(20) * t3
      t45 = ddreal(59) * t7
      t46 = s16 * t16
      t47 = ddreal(50)
      t48 = ddreal(83) * t4
      t49 = ddreal(19)
      t50 = t14 * s16
      t51 = ddreal(25) * s25
      t52 = t3 * t49
      t53 = ddreal(28) * t7
      t54 = ((t51 - t50) * s25 - t52) * s25 + t53
      t55 = t54 * s26
      t56 = t15 * s25
      t57 = t14 * t3
      t58 = ddreal(16) * t7
      t59 = ddreal(11) * t4
      t60 = t13 * t10
      t61 = ddreal(30) * s25
      t62 = ddreal(145) * s25
      t63 = ddreal(148) * s25
      t64 = ddreal(31) * s16
      t65 = ddreal(42) * t3
      t66 = (ddreal(39) * s25 + t64) * s25 + t65
      t67 = t16 + t3
      t68 = t66 * s26
      t69 = ddreal(80) * s16
      t70 = ddreal(27) * s25
      t71 = ddreal(85) * t4
      t72 = s16 * s25
      t73 = ddreal(135) * s25
      t74 = ddreal(75) * s25
      t75 = ddreal(835) * t4
      t76 = (((ddreal(588) * s16 + t74) * s25 + ddreal(1052) * t3) * s25
     & + ddreal(1150) * t7) * s25
      t77 = t76 + t75
      t78 = ddreal(654) * t7
      t79 = ((ddreal(651) * s16 + ddreal(230) * s25) * s25 + ddreal(789)
     & * t3) * s25
      t80 = t78 + t79
      t81 = (ddreal(128) * s16 + ddreal(81) * s25) * s25 + ddreal(117) *
     & t3
      t82 = ddreal(28) * s16
      t83 = -t82 - t70
      t84 = (s26)**2
      t85 = (t84)**2
      t86 = t84 * t85
      t87 = s26 * t85
      t88 = s26 * t84
      t89 = s16 * ((((ddreal(534) * s16 + t73) * s25 + ddreal(769) * t3)
     & * s25 + ddreal(798) * t7) * s25 + ddreal(525) * t4)
      t90 = ddreal(140) * t4
      t91 = ddreal(100) * s25
      t92 = ddreal(630) * t7
      t93 = ((ddreal(562) * s16 + t91) * s25 + ddreal(881) * t3) * s25
      t94 = t92 + t93
      t95 = ddreal(406) * t3
      t96 = (ddreal(470) * s16 + ddreal(200) * s25) * s25
      t97 = t95 + t96
      t98 = ddreal(106)
      t99 = t98 * s16
      t100 = ddreal(88) * s25
      t101 = t100 + t99
      t102 = s16 * (((ddreal(560) * s16 + ddreal(170) * s25) * s25 + ddr
     &eal(722) * t3) * s25 + ddreal(455) * t7)
      t103 = ddreal(78) * s25
      t104 = ddreal(120)
      t105 = ddreal(271) * t3
      t74 = (ddreal(286) * s16 + t74) * s25
      t106 = t105 + t74
      t107 = ddreal(131) * s16
      t108 = ddreal(91) * s25
      t109 = t108 + t107
      t110 = s16 * ((s25 * t104 + ddreal(299) * s16) * s25 + ddreal(231)
     & * t3)
      t111 = ddreal(57)
      t112 = t111 * s25
      t113 = ddreal(7) * s16
      t114 = t34 * s16
      t115 = s25 + t114
      t116 = ddreal(9)
      t117 = ddreal(17) * s26
      t118 = s16 * t116
      t119 = ddreal(4)
      t120 = t15 * s26
      t121 = t119 * t3
      t122 = ddreal(48) * t3
      t123 = ddreal(79) * t4
      t124 = ddreal(12)
      t125 = ddreal(10) * s25
      t126 = ddreal(95) * s16
      t127 = ddreal(36) * t3
      t128 = t124 * t7
      t129 = t15 * t54
      t130 = ddreal(15) * s25
      t131 = ddreal(176) * s16
      t132 = ddreal(581) * t4
      t133 = (((ddreal(443) * s16 + t62) * s25 + ddreal(444) * t3) * s25
     & + ddreal(415) * t7) * s25
      t134 = ddreal(82)
      t135 = ddreal(128) * t7
      t136 = ((ddreal(97) * s16 + ddreal(74) * s25) * s25 + t134 * t3) *
     & s25 + t135
      t137 = t120 * t66
      t138 = t136 * t13
      t139 = (t132 + t133) * t24
      t140 = t81 * t13
      t141 = t80 * t24
      t142 = ddreal(35) * s26
      t143 = t97 * t14 / ddreal(2)
      t144 = ddreal(20) * s25
      t145 = ddreal(40) * s16
      t146 = ddreal(21) * s16
      t147 = ddreal(10) * s26
      t148 = (s34)**2
      t149 = (t148)**2
      t150 = (t149)**2
      t151 = s34 * t148
      t152 = t151 * t149
      t153 = t148 * t149
      t154 = s34 * t149
      t155 = ((((-t147 - t113) * t148 - ddreal(76) * t88 - t110 - s26 * 
     &(t24 * t109 * s26 + t34 * t106)) * s34 + (((t13 * t101 / ddreal(2)
     & + t142) * s26 + t143) * s26 + t34 * t94) * s26 + t102) * s34 - ((
     &(-t120 * t83 + t140) * s26 + t141) * s26 + t34 * t77) * s26 - t89)
     & * s34 + s16 * (((((ddreal(287) * s16 + t112) * s25 + ddreal(329) 
     &* t3) * s25 + ddreal(293) * t7) * s25 + ddreal(377) * t4) * s25 + 
     &ddreal(357) * t10) + (((t138 + t137) * s26 + t139) * s26 + t119 * 
     &(((((t131 + t130) * s25 + ddreal(277) * t3) * s25 + ddreal(249) * 
     &t7) * s25 + ddreal(287) * t4) * s25 + ddreal(320) * t10)) * s26
      t156 = t1 * (((((-s16 * t31 - t22) * s25 - t32) * s25 - ddreal(63)
     & * t7) * s25 - t33) * s25 + t10 * t30)
      t157 = t8 * s25
      t158 = ddreal(15) * s16
      t159 = t15 * t85
      t31 = t3 * t31
      t160 = (ddreal(158) * s16 + t42) * s25
      t161 = ddreal(53) * t7
      t162 = ddreal(112) * t4
      t163 = ddreal(271) * t10
      t164 = ddreal(44) * s16
      t165 = ddreal(67) * s25
      t166 = (t165 + t164) * s25
      t167 = ddreal(88) * t7
      t168 = ddreal(139) * t4
      t169 = ((-t57 + t166) * s25 - t167) * s25 + t168
      t170 = ddreal(118)
      t171 = ddreal(123)
      t172 = ddreal(18) * t7
      t173 = ddreal(10) * t55
      t174 = ddreal(352) * s16
      t175 = ddreal(10) * t68
      t176 = -ddreal(10) * t83 * s26
      t177 = ddreal(70) * s26
      t178 = ddreal(60)
      t179 = ddreal(10) * t1 * t40
      t180 = ddreal(84) * t4
      t181 = t3 * t116
      t182 = ddreal(18) * s16
      t183 = t9 * s25
      t184 = t24 * s26
      t185 = ddreal(58)
      t186 = t7 * t185
      t187 = ddreal(42) * t7
      t188 = t149 + t4
      t189 = ddreal(23) * t7
      t190 = t2 * s26
      t191 = t18 + t7
      t192 = (t41)**2
      t77 = (((-((t120 + t113) * s26 + t121) * s34 + ((ddreal(30) * t115
     & + t117) * s26 + t118 * (t56 + t113)) * s26 + t3 * (ddreal(23) * s
     &25 + t82)) * s34 - t3 * ((ddreal(112) * s16 + t112) * s25 + ddreal
     &(84) * t3) - t49 * t85 - ((s26 * t109 + t106) * s26 + t110) * s26)
     & * s34 + t3 * (((ddreal(190) * s16 + t103) * s25 + ddreal(215) * t
     &3) * s25 + ddreal(140) * t7) + ddreal(7) * t87 + (((s26 * t101 + t
     &97) * s26 + t94) * s26 + t102) * s26) * s34 + t3 * ((((-ddreal(168
     &) * s16 - ddreal(62) * s25) * s25 - ddreal(218) * t3) * s25 - ddre
     &al(200) * t7) * s25 - t90) - (((-t83 * t84 + t80) * s26 + t77) * s
     &26 + t89) * s26 - t34 * t81 * t85
      t52 = (s25 * (((t22 - t118) * s25 + t52) * s25 - t189) + ((s34 * t
     &83 + t66) * s34 - t54) * s34 + ddreal(7) * t188) * s56 + ((t124 * 
     &t18 + t187) * s25 - t180) * s25 - t154 * t49 + ddreal(30) * t10 + 
     &((((t142 + t100 + t99) * s34 - s16 * (ddreal(234) * s16 + ddreal(1
     &40) * s26) - (ddreal(256) * s16 + ddreal(162) * s25 + ddreal(135) 
     &* s26) * s25) * s34 + ((ddreal(194) * s16 + t63) * s25 + ddreal(16
     &4) * t3) * s25 + t137 + ddreal(256) * t7) * s34 + ((t57 - t166) * 
     &s25 + t167) * s25 - t168 - t129 * s26) * s34 + t120 * t1 * t40
      t52 = t52 * s56 + (((((t22 + t26) * s25 - ddreal(28) * t3) * s25 +
     & t186) * s25 + t34 * t4) * s25 - ddreal(117) * t10) * s25 + (t179 
     &* s26 + ddreal(24) * t36) * s26 + t30 * t9 + ddreal(17) * t153 + (
     &(((-(ddreal(76) * s26 + t107 + t108) * s34 + (ddreal(424) * s16 + 
     &ddreal(352) * s25 + t177) * s26 + t95 + t96) * s34 - (t176 + t140)
     & * s26 - t78 - t79) * s34 + (t138 + t175) * s26 + t132 + t133) * s
     &34 + (((-t31 - t160) * s25 + t161) * s25 + t162) * s25 + (-t119 * 
     &t169 - t173) * s26 - t163) * s34
      t31 = t52 * s56 - t15 * t152 + ddreal(46) * t8 + ((((((s16 * t178 
     &+ ddreal(51) * s26 + t61) * s34 - (ddreal(393) * s16 + ddreal(273)
     & * s25 + ddreal(114) * s26) * s26 - t105 - t74) * s34 + ((t14 * t1
     &01 + t177) * s26 + t143) * s26 + t92 + t93) * s34 - ((t124 * t81 +
     & t176) * s26 + t141) * s26 - t75 - t76) * s34 + ((((t61 + t174) * 
     &s25 + ddreal(554) * t3) * s25 + ddreal(498) * t7) * s25 + ddreal(5
     &74) * t4) * s25 + ((t124 * t136 + t175) * s26 + t139) * s26 + ddre
     &al(640) * t10) * s34 + (((((-s16 * t170 - t56) * s25 - ddreal(121)
     & * t3) * s25 - t172) * s25 + t171 * t4) * s25 + ddreal(28) * t10) 
     &* s25 + ((-t14 * t169 - t173) * s26 - t24 * ((((t31 + t160) * s25 
     &- t161) * s25 - t162) * s25 + t163)) * s26 - ddreal(265) * t9) * s
     &34 + t179 * t88 + ddreal(36) * t36 * t84 + t180 * t18 - ddreal(56)
     & * t10 * t16 - ddreal(24) * t7 * t17 + t181 * t20 + t182 * t19 - d
     &dreal(77) * t183 + t184 * t156
      t31 = t31 * s56 + t24 * ((s16 * (t130 + t146) + (t144 + t117 + t14
     &5) * s26) * t153 + t156 * t84) + t3 * (t16 * ((((-t158 + t25) * s2
     &5 + t57) * s25 + ddreal(39) * t7) * s25 - ddreal(45) * t4) + ddrea
     &l(21) * t9) + (t155 * s34 - t41 * ((((((t125 + t126) * s25 + t127)
     & * s25 - t128) * s25 - ddreal(81) * t4) * s25 + ddreal(23) * t10) 
     &* s25 + ddreal(133) * t9 + s26 * (((((ddreal(51) * s16 + ddreal(26
     &8) * s25) * s25 + t57) * s25 - ddreal(257) * t7) * s25 + ddreal(41
     &6) * t4) * s26 + ((((ddreal(206) * s16 + ddreal(141) * s25) * s25 
     &+ t122) * s25 - ddreal(165) * t7) * s25 - t123) * s25 + ddreal(397
     &) * t10) + t129 * t88)) * s34 + ddreal(24) * t36 * t88 - ddreal(24
     &) * t157 + t159 * t1 * t40 + t114 * t1 * t29 * s26
      t52 = s16 * t1
      t29 = s34 * ((t77 * s34 + t41 * (ddreal(84) * t4 * t67 + (((((ddre
     &al(155) * s16 + t63) * s25 + ddreal(133) * t3) * s25 + ddreal(214)
     & * t7 + t68) * s26 + (((ddreal(295) * s16 + t62) * s25 + ddreal(28
     &9) * t3) * s25 + ddreal(282) * t7) * s25 + ddreal(367) * t4) * s26
     & + ((((ddreal(207) * s16 + t61) * s25 + ddreal(259) * t3) * s25 + 
     &ddreal(209) * t7) * s25 + ddreal(292) * t4) * s25 + ddreal(273) * 
     &t10) * s26 + t72 * (t16 * ((t70 + t69) * s25 + ddreal(70) * t3) + 
     &t71))) * s34 - t192 * ((((((t56 + t43) * s25 + t57) * s25 - t58) *
     & s25 - t59) * s25 + t60) * s25 + ddreal(28) * t9 + ((t14 * t46 * t
     &1 + (-t47 * t7 + ddreal(67) * t18) * s25 + t48 + t55) * s26 + ((((
     &t43 + t42) * s25 + t44) * s25 - t45) * s25 + t23) * s25 + ddreal(7
     &7) * t10) * s26)) + t31 * s56 + t52 * t29 * t84 + t1 * t40 * t87 +
     & t156 * t88
      t31 = s16 * t13
      t40 = t116 * s25
      t42 = ddreal(10) * t7
      t54 = t4 * t16
      t55 = -t34 * (s34 * t2 + t72) + t16 + t148 + t3
      t62 = s16 + s26 + s56
      t63 = t34 * s26
      t66 = (-s34 + t2 + s26) * s26 + (t63 + s16 + s25 - s34 + s56) * s5
     &6 + t72
      t68 = t24 * s25
      t74 = t119 * t4
      t75 = ddreal(7) * t10
      t76 = t15 * t7
      t77 = t4 * t49
      t78 = t24 * t52 * t18
      t79 = ((t14 * t18 + t76) * s25 - t77) * s25
      t80 = -t11 * t34 - t18 + t7
      t81 = ddreal(22)
      t83 = ddreal(42) * t10
      t89 = ((t16 * (-s25 * t81 - t43) + t128) * s25 + t14 * t4) * s25
      t92 = ddreal(32) * s25
      t93 = t7 * t49
      t94 = (-t92 - t118) * t16
      t95 = ddreal(40) * t4
      t96 = (t93 + t94) * s25 - t95
      t97 = t191 * t34 - t11
      t99 = t34 * s25
      t100 = t124 * t4
      t101 = ddreal(105) * t4
      t102 = (((ddreal(96) * s16 + ddreal(65) * s25) * s25 + ddreal(85) 
     &* t3) * s25 + ddreal(86) * t7) * s25
      t105 = t3 * t111
      t106 = ddreal(80) * t7
      t107 = ((ddreal(67) * s16 + ddreal(68) * s25) * s25 + t105) * s25
      t108 = t106 + t107
      t109 = t14 * t67 + t72 * t15
      t110 = t124 * s25
      t117 = t4 * t47
      t120 = t24 * t88
      t129 = ddreal(70) * t7
      t132 = ((s25 * t47 + ddreal(87) * s16) * s25 + ddreal(92) * t3) * 
     &s25
      t133 = ddreal(72) * s25
      t136 = ddreal(80) * t3
      t137 = (t126 + t133) * s25
      t138 = t136 + t137
      t139 = ddreal(110)
      t140 = ddreal(86) * s16
      t141 = ddreal(70) * t4
      t142 = ddreal(105) * t3
      t143 = (ddreal(141) * s16 + ddreal(85) * s25) * s25
      t155 = ddreal(38) * s25
      t156 = ddreal(115)
      t160 = ddreal(40) * s25
      t161 = t3 * t156
      t163 = s16 * t41
      t166 = ddreal(59) * s16
      t167 = ddreal(14) * s16
      t168 = t16 * t67
      t169 = ddreal(11) * t3
      t173 = t24 * t16
      t175 = t13 * t7
      t176 = ddreal(18) * s26
      t177 = t176 * t97
      t179 = t116 * s26
      t180 = t179 * t109
      t193 = t49 * s25
      t194 = t13 * s26
      t195 = t20 * t2
      t196 = t124 * s16
      t197 = t3 * t13
      t198 = t10 - t154
      t199 = t18 * t1
      t200 = t119 * t2
      t93 = (t13 * t198 + t78 + t79 - (((-(t179 + t145 + t155) * s34 + t
     &136 + t137 + ddreal(36) * t190) * s34 - t106 - t107 - t180) * s34 
     &+ (-t93 - t94) * s25 + t177 + t95) * s34 + t179 * t1 * t80) * s56 
     &+ t195 * t24 - t34 * (t151 * (t148 * (t193 + t194 + t146) + (ddrea
     &l(18) * t190 + t138) * s26 + t129 + t132) - t1 * (s16 * ((-t173 - 
     &t169) * s25 + t175) - t14 * t168) * s26) + ((t148 * ((ddreal(76) *
     & s25 + t69 + t179) * s26 + t142 + t143) + (t108 * t34 + t180) * s2
     &6 + t101 + t102) * s34 + (t34 * t96 - t177) * s26 - t83 + t89) * s
     &34 + ddreal(7) * t9 + ddreal(7) * t153 + t116 * t1 * t80 * t84 + t
     &11 * (((-t22 + t196) * s25 - t197) * s25 - ddreal(11) * t7)
      t94 = (s56)**2
      t106 = s56 * t94
      t80 = (-t200 * t151 + t11 * (-t37 + t99) - ((-t148 - t109) * s34 +
     & t34 * t97) * s34 + t4 - t199) * t106 + t1 * t80 * t88 - t52 * t20
      t60 = t80 * t24 + t34 * (t8 - t152) + ((t78 + t79 + t60) * s26 + t
     &1 * (((t16 * (-t68 - t50) - t128) * s25 - t74) * s25 + t75)) * s26
     & + ((((((ddreal(7) * s26 + t167 + t110) * s34 - (ddreal(38) * s26 
     &+ t61 + t166) * s25 - t13 * t84 - ddreal(42) * t163) * s34 + ((ddr
     &eal(99) * s16 + t160) * s25 + t161) * s25 + t120 + t129 + s26 * ((
     &t145 + t155) * s26 + t142 + t143)) * s34 - ((t190 * t124 + t138) *
     & s26 + t34 * (t132 + t129)) * s26 - t12 * t139 - t18 * (t61 + t140
     &) - t141) * s34 + ((((t110 + t164) * s25 + ddreal(40) * t3) * s25 
     &+ ddreal(36) * t7) * s25 + t117) * s25 + t83 + s26 * (s26 * t108 +
     & t101 + t102) + t120 * t109) * s34 + ((t18 * (-t99 - t158) + t100)
     & * s25 - t75) * s25 - ddreal(14) * t9 + s26 * (s26 * t96 - t83 + t
     &89) - t97 * t14 * t88) * s34 + t93 * s56 - t183 + t54 * (t22 - t11
     &3)
      t78 = -t116 * s25 * (t16 * (s16 + s25) + t7) + t23
      t79 = ddreal(45) * s16
      t80 = (t16 * (-t40 - t79) - t189) * s25 + ddreal(17) * t4
      t89 = t3 * t34
      t93 = (-t158 - t25) * s25
      t96 = s16 - t68
      t97 = ddreal(53) * s25
      t101 = ddreal(20) * t7
      t99 = s16 + t99
      t102 = ddreal(11) * t16
      t107 = t88 * t99
      t108 = ddreal(147)
      t109 = ddreal(30) * t3
      t113 = t119 * s25
      t120 = ddreal(340) * t7
      t132 = ((s25 * t170 + ddreal(314) * s16) * s25 + ddreal(357) * t3)
     & * s25
      t136 = t120 + t132
      t137 = ddreal(91) * s16
      t138 = ddreal(63) * s25
      t142 = (t138 + t137) * s25 + ddreal(76) * t3
      t143 = ddreal(20) * s16 + t70
      t155 = ddreal(245) * t7
      t177 = ((ddreal(202) * s16 + t160) * s25 + ddreal(301) * t3) * s25
      t180 = ddreal(235) * t3
      t201 = (ddreal(261) * s16 + ddreal(122) * s25) * s25
      t202 = t180 + t201
      t203 = ddreal(73) * s16
      t204 = t203 + t165
      t205 = ddreal(88) * s16
      t206 = ddreal(28) * s25
      t207 = ddreal(35) * s16
      t208 = ddreal(40) * t7
      t209 = ddreal(26) * s16
      t210 = ddreal(11) * s25
      t211 = t7 * t14
      t212 = ddreal(13) * t4
      t213 = ddreal(47) * s16
      t214 = (t213 + t112) * s25 + t109
      t215 = (((ddreal(465) * s16 + ddreal(348) * s25) * s25 + ddreal(36
     &9) * t3) * s25 + ddreal(474) * t7) * s26
      t216 = ddreal(15) * t85
      t217 = ddreal(152) * t7
      t218 = t24 * t84
      t219 = ddreal(20) * s26
      t45 = (t148 * (t148 * (ddreal(63) * s25 * s26 + t140 * s26 + ddrea
     &l(20) * t16 + ddreal(63) * t3 + ddreal(65) * t72 + ddreal(21) * t8
     &4) + (((t125 + t205) * s25 + ddreal(143) * t3) * s25 + t217) * s25
     & + t90 + s26 * (t84 * t34 * t143 + t136) + t218 * t142) + t41 * ((
     &(((t209 + t210) * s25 + t197) * s25 - t211) * s25 - t212 + ddreal(
     &15) * t107) * s25 + ddreal(35) * t10 + s26 * ((((s25 * t98 - t50) 
     &* s25 + t169) * s25 + t208) * s26 + (((ddreal(17) * s16 + t103) * 
     &s25 + t122) * s25 - t45) * s25 + t48))) * s34
      t48 = t1 * ((t93 - t89) * s25 + t76)
      t76 = ddreal(36) * s16
      t197 = ddreal(24) * t3
      t220 = t14 * s26
      t221 = s16 * t191
      t222 = ddreal(14) * t3
      t223 = ddreal(26) * t16
      t224 = ddreal(41) * t221
      t225 = ddreal(13) * t3
      t226 = ((t196 + t97) * s25 + t225) * s25 + t101
      t227 = ddreal(25) * t3
      t228 = t170 * t10
      t229 = t61 * t88
      t230 = t14 * t84
      t231 = t3 * t16
      t232 = ddreal(201) * s25
      t233 = ddreal(30) * s26
      t234 = ddreal(48) * s25
      t235 = t184 * t1
      t236 = t58 * s25
      t237 = t61 * t84
      t238 = t119 * s26
      t239 = ddreal(40) * s26
      t240 = t178 * s26
      t241 = ddreal(13) * s16
      t242 = ddreal(15) * s26
      t65 = (((s16 * (t167 + t110) + (ddreal(13) * s26 + t144 + t207) * 
     &s26) * s34 - s16 * ((t206 + t166) * s25 + t65) - ((ddreal(14) * s2
     &6 + t138 + t140) * s26 + t34 * ((ddreal(65) * s16 + t144) * s25 + 
     &ddreal(63) * t3)) * s26) * s34 + s16 * (((t92 + t205) * s25 + t161
     &) * s25 + t129) + t159 + ((s26 * t204 + t202) * s26 + t155 + t177)
     & * s26) * s34 - s16 * ((((ddreal(66) * s16 + t25) * s25 + ddreal(8
     &6) * t3) * s25 + t139 * t7) * s25 + t141) - ((t143 * t84 + t136) *
     & s26 + (((t144 + t131) * s25 + ddreal(286) * t3) * s25 + ddreal(30
     &4) * t7) * s25 + ddreal(280) * t4) * s26 - t34 * t142 * t88
      t61 = (s25 * (((-t37 + t25) * s25 - t225) * s25 - t35 + t242 * t96
     & * t1) + t15 * t188 - ((t143 * s34 - (t242 + t213 + t112) * s25 - 
     &t109) * s34 + (s16 * (t233 + t241) + (t196 + t240 + t97) * s25) * 
     &s25 + t101) * s34) * s56 + ((((t40 + t76) * s25 - ddreal(45) * t3)
     & * s25 + t189) * s25 - t95) * s25 - t34 * s34 * (t34 * t226 * s26 
     &+ t61 * t99 * t84 + t148 * (s16 * (ddreal(76) * s16 + t239) + (ddr
     &eal(54) * s26 + t137 + t138) * s25) + ddreal(26) * t17 + t224 + dd
     &real(14) * t231 - t236) + ddreal(17) * t10 - ddreal(14) * t154 + t
     &148 * (t11 * t171 + t148 * (t219 + t203 + t165) + t238 * t214 + dd
     &real(116) * t18 + t237 + ddreal(155) * t46 + ddreal(158) * t7) + t
     &237 * t96 * t1 + t238 * t48
      t58 = t61 * s56 + t10 * (-t234 + t146) + ddreal(13) * t153 + ddrea
     &l(27) * t46 * ((t16 - t3) * s25 + t7) - ((-ddreal(216) * s16 * t18
     & + (((ddreal(42) * s26 + t140 + t138) * s34 - (ddreal(219) * s16 +
     & t232 + t233) * s26 - t180 - t201) * s34 + t14 * s26 * (s26 * t143
     & + t142) + t120 + t132) * s34 - t111 * t17 - t230 * t214 - t73 * t
     &7 - t215 - t229 - ddreal(207) * t231 - ddreal(275) * t4) * s34 + (
     &(((ddreal(104) * s16 + t210) * s25 + t227) * s25 + t187) * s25 - d
     &dreal(72) * t4 + t107 * t178) * s25 + t14 * s26 * (t226 * s26 + s2
     &5 * ((t222 + t223) * s25 - t58) + t224) + t228) * s34 + t230 * t48
     & + t229 * t96 * t1 + t235 * t80
      t45 = t58 * s56 - t119 * (-t48 * t88 + t152) + t3 * (((((t70 - t76
     &) * s25 + t197) * s25 - t35) * s25 - t77) * s25 + ddreal(11) * t10
     &) - t34 * t45 + (t148 * (t148 * (ddreal(26) * s26 + t144 + t207) +
     & ((t204 * t24 + t219) * s26 + t34 * t202) * s26 + t155 + t177) + (
     &((ddreal(136) * t3 + ddreal(91) * t72) * s25 + ddreal(119) * t7) *
     & s25 + t4 * t98 + t216) * s25 + ddreal(189) * t10 + t119 * (t214 *
     & t88 + t20) + s26 * ((((ddreal(432) * s16 + ddreal(114) * s25) * s
     &25 + ddreal(414) * t3) * s25 + ddreal(270) * t7) * s25 + ddreal(55
     &0) * t4 + t215)) * t148 + t218 * t1 * t80 + t220 * t52 * t78 + t13
     &0 * t96 * t1 * t85
      t58 = t96 * t1
      t61 = t58 * s25 * t87 + t52 * t78 * t84 + s25 * (-t119 * s25 * (s1
     &6 + s34) + t173 + t148 + t3 - t114 * s34) * s56 * (t94)**2
      t73 = t7 * t16
      t77 = s16 + s26 - s34 + s56
      t78 = ddreal(11) * s16
      t95 = ddreal(14) * s25
      t98 = (t95 - t78) * s25 + t27
      t120 = ddreal(90)
      t130 = ddreal(64) * t3
      t131 = t16 - t3
      t132 = -t14 * s16 * t131 + (t16 * t34 - t181) * s25
      t136 = ddreal(38) * s16
      t137 = ddreal(77) * t7
      t138 = t7 * t1
      t140 = ddreal(84) * s25
      t141 = ddreal(10) * s16
      t39 = (((-t95 - t126) * s25 + t39) * s25 - t7) * s25 + t74
      t126 = t132 * s26
      t142 = ddreal(18) * t3
      t143 = (t50 + t210) * s25 + t142
      t161 = ddreal(38) * t3
      t165 = t120 * t7
      t166 = ((t234 - t207) * s25 - t161) * s25
      t167 = t47 * t3
      t177 = ddreal(55) * s16
      t180 = ddreal(55) * t3
      t181 = ddreal(28) * t16
      t187 = ddreal(34) * t7
      t189 = t56 + t50
      t201 = s16 * (((-t181 + t180) * s25 - t187) * s25 + t23)
      t202 = t24 * t7
      t204 = (((-ddreal(39) * s16 + t22) * s25 + t222) * s25 - ddreal(47
     &) * t7) * s25 + ddreal(20) * t4
      t213 = t132 * t88
      t214 = t10 * s25
      t215 = -ddreal(264) * s25
      t219 = (((ddreal(460) * s16 + t215) * s25 + ddreal(492) * t3) * s2
     &5 + ddreal(608) * t7) * s25 - ddreal(232) * t4
      t222 = t166 + t165
      t224 = ddreal(10) * t222
      t225 = s26 * t143
      t226 = ddreal(278) * t10
      t229 = ((((t158 - t92) * s25 + ddreal(372) * t3) * s25 + ddreal(33
     &5) * t7) * s25 + ddreal(54) * t4) * s25
      t237 = ddreal(76) * t7
      t243 = ((t146 + t133) * s25 + t169) * s25
      t244 = t237 + t243
      t245 = (ddreal(36) * s25 + t207) * s25 + t167
      t246 = ddreal(20) * t245
      t247 = t244 * t13
      t248 = ddreal(168) * t3
      t249 = (ddreal(179) * s16 + ddreal(156) * s25) * s25
      t250 = t248 + t249
      t251 = t177 + t234
      t252 = -ddreal(133) * s16 - ddreal(160) * s25
      t253 = t250 * t119
      t254 = (ddreal(153) * s16 + ddreal(140) * s25) * s25
      t255 = ddreal(36) * s26
      t256 = t98 * t88
      t257 = ddreal(10) * t201
      t258 = ddreal(74) * t4
      t259 = ddreal(30) * t126
      t260 = ddreal(66) * s25
      t242 = t242 * t143
      t261 = ddreal(240) * s25
      t262 = ddreal(132)
      t263 = ddreal(35) * s25
      t264 = ddreal(45) * s25
      t265 = (((-s16 * t134 - t264) * s25 + ddreal(265) * t3) * s25 + dd
     &real(490) * t7) * s25 + ddreal(670) * t4
      t266 = (s25 * t252 + ddreal(10) * t3) * s25 + t129
      t262 = s16 * ((((s16 * t262 - t263) * s25 + ddreal(423) * t3) * s2
     &5 + ddreal(570) * t7) * s25 + ddreal(630) * t4)
      t267 = ddreal(78) * s16
      t268 = ddreal(180)
      t268 = (s16 * t268 - t125) * s25
      t269 = ddreal(419) * t3
      t270 = ddreal(448) * t7
      t271 = ddreal(48) * s16
      t272 = t178 * s25
      t273 = ((t272 + t271) * s25 - ddreal(185) * t3) * s25
      t274 = ddreal(324) * t7
      t275 = -t274 + t273
      t276 = t130 + t254
      t277 = t205 + t140
      t47 = (((-s16 * t47 - t206) * s25 + ddreal(32) * t3) * s25 - ddrea
     &l(27) * t7) * s25 + t212
      t203 = (t203 + t125) * s25 + ddreal(98) * t3
      t49 = (-s16 * t49 - t264) * s25
      t212 = ddreal(68) * t3
      t264 = t212 + t49
      t278 = -ddreal(54) * s16 - ddreal(64) * s25
      t279 = ddreal(56) * t4
      t280 = t119 * s16
      t281 = ddreal(211)
      t282 = ddreal(30) * s16
      t283 = ddreal(7) * s25
      t284 = ddreal(438)
      t285 = ddreal(96) * t10
      t286 = ddreal(382) * t7
      t287 = ddreal(241)
      t288 = (-t25 - t205) * s25
      t289 = ddreal(325) * t3
      t290 = ddreal(558) * t7
      t291 = ddreal(360) * t4
      t292 = ddreal(772) * t10
      t293 = (((ddreal(63) * s16 + ddreal(300) * s25) * s25 - ddreal(738
     &) * t3) * s25 - ddreal(684) * t7) * s25 - ddreal(780) * t4
      t294 = t24 * t266
      t295 = t24 * t276
      t296 = s25 + s26
      t51 = ((t148 * (s34 * (t21 - t220) + (-t193 + t194 - t271) * s16 +
     & ddreal(36) * s26 * t296) + s16 * ((-t269 - t268) * s25 - t270) + 
     &t178 * t85 + ((t119 * t277 * s26 + t295) * s26 + t34 * t275) * s26
     &) * s34 + ((((-t15 * t251 - t176) * s26 - t253) * s26 + t294) * s2
     &6 + t265 * t34) * s26 + t262) * s34 + s16 * (((((-ddreal(56) * s16
     & + t51) * s25 - ddreal(299) * t3) * s25 - t286) * s25 - ddreal(345
     &) * t4) * s25 - ddreal(560) * t10) + (-t34 * ((((t289 + t288) * s2
     &5 + t290) * s25 + t291) * s25 + t292) + t84 * (t147 * t245 + t247)
     &) * s26 + t84 * (t88 * t124 * t189 + t293)
      t28 = s34 * ((t51 * s34 + t41 * (-t14 * (t143 * t85 + t195) + t3 *
     & ((((ddreal(184) * s16 + t232) * s25 + ddreal(108) * t3) * s25 + d
     &dreal(14) * t7) * s25 + ddreal(308) * t4) + (((((s16 * t287 - t261
     &) * s25 + ddreal(226) * t3) * s25 - ddreal(342) * t7) * s26 + (((d
     &dreal(700) * s16 + t215) * s25 + ddreal(251) * t3) * s25 + t286) *
     & s25 + t139 * t4) * s26 + ((((ddreal(309) * s16 - ddreal(96) * s25
     &) * s25 + ddreal(416) * t3) * s25 + ddreal(754) * t7) * s25 - ddre
     &al(220) * t4) * s25 + ddreal(724) * t10) * s26)) * s34 - t192 * (s
     &16 * (((((-s16 * t30 + ddreal(102) * s25) * s25 + ddreal(87) * t3)
     & * s25 - t172) * s25 - t123) * s25 + t285) + (((((s16 * t284 - t27
     &2) * s25 - ddreal(284) * t3) * s25 + ddreal(254) * t7) * s25 - t27
     &9) * s26 + t119 * (((((s16 * t139 - t110) * s25 - ddreal(99) * t3)
     & * s25 + t170 * t7) * s25 - ddreal(80) * t4) * s25 + t28)) * s26 -
     & t124 * (t19 + t213)))
      t51 = t138 * ((((-t141 - t140) * s25 + ddreal(26) * t3) * s25 - dd
     &real(51) * t7) * s25 + ddreal(29) * t4)
      t123 = s16 * t47
      t139 = ddreal(26) * s25
      t170 = ddreal(600)
      t49 = (((((t280 + t25 + t255) * s34 + (-ddreal(162) * s16 - ddreal
     &(192) * s25 - ddreal(108) * s26) * s26 + t212 + t49) * s34 + t104 
     &* t88 + t273 - t274 + s26 * (t14 * t277 * s26 + t295)) * s34 + s16
     & * (((-s25 * t134 + ddreal(265) * s16) * s25 + ddreal(490) * t3) *
     & s25 + ddreal(670) * t7) + ((-t14 * t250 - t147 * t251) * s26 + t2
     &94) * s26 - ddreal(45) * t17 - ddreal(45) * t85) * s34 + (((-t289 
     &- t288) * s25 - t290) * s25 - t291) * s25 + (((t233 * t189 + t246)
     & * s26 + t124 * t244) * s26 + t293) * s26 - t292) * s34 + s16 * ((
     &(((ddreal(255) * s16 - ddreal(51) * s25) * s25 + ddreal(300) * t3)
     & * s25 + ddreal(431) * t7) * s25 - ddreal(103) * t4) * s25 + ddrea
     &l(516) * t10) + t24 * ((t226 + t229) * s26 - t19) + t84 * ((-t224 
     &- t242) * s26 + ddreal(3) / ddreal(2) * t219)
      t82 = (-t143 * t148 - t188 * t24 + t34 * s34 * (t148 * t189 + t132
     &) + ddreal(14) * t72 * t67 - ddreal(25) * t231) * s56 + s16 * (((-
     &t180 + t181) * s25 + t187) * s25 - t23 - t220 * t1 * t98) + t151 *
     & (t124 * t148 + t34 * (ddreal(36) * s16 * s26 + (t233 + t207) * s2
     &5 + ddreal(36) * t16 + t167)) - ((t220 * t143 + t148 * (t176 + t17
     &7 + t234) + t165 + t166) * s34 - s16 * (((t82 - t103) * s25 - ddre
     &al(94) * t3) * s25 + t208) - t124 * (t17 + t126)) * s34
      t69 = t82 * s56 + s16 * (t10 * t119 + t16 * (((ddreal(81) * s16 + 
     &t95) * s25 - ddreal(111) * t3) * s25 + ddreal(17) * t7) - ddreal(1
     &5) * t1 * t98 * t84) - t15 * (s26 * t201 + t214) + t151 * (t34 * (
     &(ddreal(15) * t189 * s26 + t15 * t245) * s26 + t237 + t243) - ddre
     &al(18) * t151) + ((t148 * ((t240 + t205 + t140) * s34 - (ddreal(27
     &5) * s16 + ddreal(45) * s26 + t261) * s26 - t248 - t249) + (((s16 
     &* t156 - t260) * s25 + t171 * t3) * s25 + t217) * s25 + (-t15 * t2
     &22 - t242) * s26 - t185 * t4) * s34 + s25 * ((((t110 - t69) * s25 
     &- ddreal(114) * t3) * s25 + t211) * s25 - t258) + (ddreal(10) * t2
     &04 + t259) * s26) * s34
      t69 = t69 * s56 + t153 * (s34 * t124 - t34 * (t255 + t92 + t26)) +
     & ddreal(26) * t8 + ((((((s26 * t104 + ddreal(336) * s25 + t174) * 
     &s26 + t130 + t254) * s34 + (-t178 * t84 - t253) * s26 + t16 * t252
     & - ddreal(10) * t251 * t84 + ddreal(10) * t11 + t129) * s34 + (((t
     &146 + t91) * s25 - ddreal(246) * t3) * s25 - ddreal(228) * t7) * s
     &25 + ((t239 * t189 + t246) * s26 + t247) * s26 - ddreal(260) * t4)
     & * s34 + ((-ddreal(20) * t225 - t224) * s26 + t219) * s26 + t226 +
     & t229) * s34 - t114 * t20 + t119 * t19 - ddreal(256) * t3 * t17 + 
     &ddreal(89) * t7 * t18 + ddreal(20) * t204 * t84 + ddreal(40) * t21
     &3 + ddreal(155) * t214 - ddreal(214) * t54 - ddreal(136) * t9 + t1
     &94 * s25 * ((((t22 - t145) * s25 - t105) * s25 + t202) * s25 - ddr
     &eal(37) * t4)) * s34 + t238 * t52 * t39 - ddreal(20) * t256 * t52 
     &- t257 * t84 + t228 * t16 - ddreal(164) * t4 * t18 + ddreal(56) * 
     &t3 * t20 + ddreal(44) * t7 * t17 - ddreal(80) * t183
      t8 = (t69 * s56 + t14 * t52 * s26 * (s26 * t39 + t123) - t150 * t2
     &4 + ddreal(29) * t5 + (t49 * s34 - t41 * ((((((t267 - t110) * s25 
     &+ t248) * s25 - ddreal(111) * t7) * s25 + ddreal(218) * t4) * s25 
     &- ddreal(239) * t10) * s25 + ((-ddreal(10) * (((-ddreal(84) * s16 
     &+ t110) * s25 + ddreal(46) * t3) * s25 - ddreal(67) * t7) * s25 - 
     &ddreal(10) * t4 * t81 - t259) * s26 + ((((s16 * t170 - t133) * s25
     & - ddreal(156) * t3) * s25 + ddreal(424) * t7) * s25 - ddreal(226)
     & * t4) * s25 + ddreal(220) * t10) * s26 + ddreal(188) * t9)) * s34
     & - t216 * t52 * t98 - t257 * t88 + ddreal(84) * t7 * t20 - t258 * 
     &t17 + ddreal(77) * t9 * t16 - ddreal(36) * t10 * t18 - ddreal(80) 
     &* t157) * s56 + t14 * t52 * t84 * (t123 - t256) + t34 * (t153 * (s
     &16 * t203 - ddreal(36) * t88 + s26 * (ddreal(3) / ddreal(2) * s26 
     &* t278 + t264)) + t51 * s26) + t28 + ddreal(13) * t6 + t52 * t119 
     &* t39 * t88 - t159 * t201 - t285 * t17 + t178 * t9 * t18 + t279 * 
     &t20 - ddreal(7) * t8 * t16 - t139 * t5
      t28 = t41 * t192 * ((t16 * (t16 * (-t113 + t136) + t137) + t75) * 
     &s25 + t34 * s26 * (-t84 * t132 + ((((s16 * t185 - t22) * s25 - ddr
     &eal(72) * t3) * s25 + ddreal(75) * t7) * s25 - t117) * s25 + ddrea
     &l(24) * t10) + t84 * ((((s16 * t120 - t110) * s25 - t130) * s25 + 
     &t208) * s25 - t74) + ddreal(16) * t9 - ddreal(56) * t231 * t67) * 
     &s34
      t33 = t192 * ((((((-t68 + t241) * s25 + ddreal(27) * t3) * s25 + d
     &dreal(93) * t7) * s25 - t33) * s25 + ddreal(63) * t10) * s25 + ddr
     &eal(56) * t9 + (((((s16 * t111 - t234) * s25 + t167) * s25 - ddrea
     &l(54) * t7 - t225) * s26 + (((s16 * t281 - t260) * s25 + t44) * s2
     &5 + t186) * s25 + ddreal(68) * t4) * s26 + ((((s16 * t108 - t92) *
     & s25 - t89) * s25 + ddreal(238) * t7) * s25 - t162) * s25 + ddreal
     &(196) * t10) * s26) * t148
      t49 = t41 * (s16 * (((((-t283 + t164) * s25 + ddreal(97) * t3) * s
     &25 + t208) * s25 + ddreal(175) * t4) * s25 + ddreal(112) * t10) + 
     &t34 * s26 * (-t189 * t85 + ((((-t40 + t50) * s25 + ddreal(101) * t
     &3) * s25 + ddreal(171) * t7) * s25 + t71) * s25 + ddreal(224) * t1
     &0) + t84 * ((((-ddreal(144) * s25 + t282) * s25 + t161) * s25 - dd
     &real(64) * t7) * s26 + (((s16 * t171 - t91) * s25 + ddreal(216) * 
     &t3) * s25 + ddreal(190) * t7) * s25 + ddreal(324) * t4) - t119 * (
     &t3 * t81 - t93) * t85) * t151
      t47 = t3 * (t1 * t47 * t88 + t5)
      t64 = t10 * t16 * (((t95 - t64) * s25 + ddreal(29) * t3) * s25 - d
     &dreal(13) * t7)
      t6 = t34 * t47 - s25 * t6 + t8 * s56 + t149 * (-t24 * t86 + t3 * (
     &(((t267 + t113) * s25 + ddreal(145) * t3) * s25 + t155) * s25 + t9
     &0) + ((((-t251 * s26 - t250) * s26 + t266) * s26 + t265) * s26 + t
     &262) * s26) + t152 * (t124 * t88 + t3 * (-ddreal(16) * s16 - ddrea
     &l(13) * s25) + s26 * (-s16 * (t193 + t271) + s26 * (t280 + t25))) 
     &+ t153 * (t84 * (s26 * t278 + t264) - ddreal(18) * t85 + t279 + dd
     &real(29) * t231 + t137 * s25 + t114 * t203 * s26) + t154 * (t124 *
     & t87 - ddreal(112) * t7 * t67 + (-s16 * ((t269 + t268) * s25 + t27
     &0) + ((s26 * t277 + t276) * s26 + t275) * s26) * s26 - t11 * (ddre
     &al(189) * t3 + t223)) - t87 * t201 + t4 * t1 * (t16 * ((-ddreal(56
     &) * s25 + t145) * s25 - t44) + ddreal(13) * t138) * s26 + t51 * t8
     &4 + t52 * t39 * t85 - t52 * t98 * t86 - t28 + t33 - t49 + (-s26 + 
     &t114) * (s16 + t184) * t150 + t64
      t8 = t24 * t4
      t25 = (-t210 + t31) * t18
      t28 = -t119 * t16 * t1 + t3 * (s16 - t68)
      t33 = t28 * t84
      t39 = -t113 + t21
      t47 = ddreal(15) * t221
      t35 = s25 * ((ddreal(25) * t16 + t121) * s25 + t35)
      t49 = ((ddreal(24) * s25 + t21) * s25 + t89) * s25
      t51 = t49 + t128
      t27 = (s16 + t22) * s25 + t27
      t64 = ((t79 + t160) * s25 + ddreal(39) * t3) * s25
      t68 = ddreal(30) * t7
      t69 = s16 * t81
      t71 = (t69 + t139) * s25 + t142
      t37 = t37 + t113
      t74 = t17 + t4
      t75 = t283 + t50
      t81 = ddreal(16) * s25
      t82 = t14 * t12
      t86 = t184 * t27
      t87 = (t1)**2
      t90 = t84 * t87
      t91 = s26 * t87
      t92 = t7 * s25
      t93 = t24 * (t16 * (-t231 + t90) + t153) * s25 + t34 * t18 * (t91 
     &* t99 + t92) + ((s25 * ((t34 * t51 + t86) * s26 + t35 + t47) + (-s
     &25 * (s26 * (t24 * t37 * s26 + t34 * t71) + t68 + t64) + (s25 * ((
     &t263 + t79) * s25 + (t184 + t206 + t43) * s26 + t109) - s25 * (t15
     &8 + t81 + t220) * s34) * s34) * s34) * s34 + (t73 * t14 - t63 * (-
     &t25 - t82 + t8)) * s25 - t13 * t17 * t131 - t24 * s25 * (t4 * (s16
     & - s25) + t33)) * s34 + t18 * t17
      t98 = s25 * (-t114 * t18 + (((-t37 + s34) * s34 + t27) * s34 - t28
     &) * s34 + t168) * t106
      t76 = s25 * (((t34 * t75 + s26) * s26 + t15 * ((t283 + t118) * s25
     & + t57)) * s26 + ddreal(20) * t191 + t72 * (ddreal(37) * s25 + t76
     &)) * t149
      t37 = s25 * (((s26 * t37 + t71) * s26 + t64 + t68) * s26 + ddreal(
     &15) * t74 + t72 * ((t209 + t206) * s25 + t197)) * t151
      t64 = t67 * s16
      t68 = -t89 * t20 + t153 * (t2 * t34 + s26) * t24
      t71 = s16 + s25 + s26 - s34 + s56
      t89 = t17 * t2
      t103 = t46 * t2
      t104 = -t175 * s25 + t24 * t74 + t103
      t105 = ddreal(15) * t4
      t106 = (t16 * (t40 + t31) - t202) * s25
      t111 = -t280 + t56
      t120 = t72 * t111 + t124 * t191
      t57 = -t13 * t16 - t57
      t121 = ((ddreal(53) * s16 + ddreal(42) * s25) * s25 + t122) * s25
      t122 = t178 * t7
      t123 = t116 * t67 + t72 * t13
      t58 = t58 * t115 * t84
      t67 = (((-ddreal(7) * t3 + t173) * s25 + t34 * (-t151 + t64) + s34
     & * (s34 * t75 + t57)) * s56 - t14 * (t151 * (s26 + t200 - s34) - t
     &17 - t4) + t148 * ((ddreal(21) * s25 + t182) * s26 + ddreal(36) * 
     &t67 + ddreal(32) * t72) + t34 * (-(s16 * (t111 * s25 + t118 * s26)
     & + t124 * (t16 * t296 + t7)) * s34 + t103) - t236 + t235 * t96 * t
     &115) * s56 + t119 * (s26 * t104 + t18 * t3) + t14 * t198 + t148 * 
     &((t123 * t13 + t184 * t75) * s26 + t121 + t122) + t24 * (t149 * (t
     &141 + t40 + t238) + t89 + t58) - t34 * s34 * (t148 * ((t282 + t263
     &) * s16 + ddreal(24) * t16 + ddreal(24) * t190 + t218) + t106 + s2
     &6 * (-ddreal(3) / ddreal(2) * t57 * s26 + t34 * t120) + t105) - t9
     &2 * (t56 + t78)
      t78 = ddreal(1) / s16
      t55 = ddreal(1) / (t55)**2
      t66 = ddreal(1) / (t66)**2
      t62 = ddreal(1) / (t62)**2
      t6 = t6 * MB1011(1) + t119 * t9 * t191 + t14 * (t19 * t7 + t36 * t
     &85) - t24 * t3 * (t183 - t190 * t1 * (-t12 * t13 + t18 * (t21 - t2
     &2) + t23)) + t29 + t54 * (t16 * (-t40 + t31) - t42)
      t5 = -(t24 * t61 - t34 * (t152 * (s16 + t63) - t5) + s34 * ((t65 *
     & s34 + t41 * (((((ddreal(34) * s16 + t113) * s25 + t127) * s25 + t
     &42) * s25 + t24 * t85 + t117) * s25 + t83 + ((((t112 + t164) * s25
     & + t109) * s26 + ((ddreal(98) * s16 + ddreal(116) * s25) * s25 + d
     &dreal(79) * t3) * s25 + t135) * s26 + (((ddreal(100) * s16 + t112)
     & * s25 + ddreal(109) * t3) * s25 + ddreal(56) * t7) * s25 + t108 *
     & t4) * s26)) * s34 - t192 * ((((t102 + t44) * s25 - ddreal(26) * t
     &7) * s25 + t23 + t107 * t14) * s25 + ddreal(14) * t10 + s26 * ((((
     &-t196 + t97) * s25 + t3) * s25 + t101) * s26 - ddreal(40) * t11 * 
     &t1 + t18 * (s25 * t30 - t43) + ddreal(42) * t4))) + t45 * s56 - t1
     &57 + t3 * t1 * ((((-t70 + t118) * s25 - ddreal(15) * t3) * s25 - t
     &175) * s25 + t59) * s26 + t48 * t85 + t1 * t80 * t88 + ddreal(18) 
     &* t199 * t4 + t73 * (t116 * t18 - t42)) * MB1001(0) + t6 * t78
      t6 = (t71)**2
      t7 = (t68 * s25 - s25 * t152 + (s16 * t18 * (t16 * t39 + t202) + (
     &(-t2 * (-t13 * t199 + t38 * t7) - t33) * s26 - t19 - t9) * s25 + t
     &84 * (t14 * t231 * t2 + s25 * (t25 - t8))) * s34 + t93 * s56 + t14
     &8 * (t14 * s25 * ((t17 + t4) * s25 + t10) + s25 * (((s26 * t27 + t
     &51) * s26 + t35 + t47) * s26 + t73) + s16 * t17 * (t114 + t110)) +
     & t94 * (t19 * t34 - t24 * (t16 * (-t91 + t46) + t154) * s25 + s25 
     &* (((((t184 + t95 + t196) * s34 - s16 * (t179 + t182) - (s26 * t12
     &4 + t139 + t69) * s25) * s34 + t128 + t49 + t86) * s34 - t24 * (s2
     &6 * t28 + t4) + t25 + t82) * s34 + t73)) + t91 * t17 * t115 + t90 
     &* t18 * t99 + t98 - s25 * ((t184 + t81) * s26 + ddreal(15) * t16 +
     & ddreal(15) * t163 + ddreal(24) * t72) * t154 + t76 - t37 + t87 * 
     &t18 * t88 + t64 * t20) * MB1101(0) * t55
      t6 = t78 * (t66 * (-t77 * t71 * t6 * MB1111D2(0) + t7) + (t77)**2 
     &* MB1011(0) * t62)
      result = -t66 * t55 * (t60 * I300s16(1) + t5 * t62 - s25 * t71 * (
     &-t15 * t73 * t1 - t24 * t52 * t17 - t34 * (-t1 * t10 - t104 * t84 
     &+ t154 * (t184 + t56 + t50 - s34)) + (t58 + t24 * t89 + t3 * ((-s2
     &5 * t39 - t169) * s25 + t211)) * s26 + (((((ddreal(27) * s26 + t14
     &4 + t136) * s25 + ddreal(30) * t163 + t230) * s34 - ((t144 + t79) 
     &* s25 + t3 * t30) * s25 - ((ddreal(24) * t2 + t63) * s26 + (ddreal
     &(70) * s16 + t234) * s25 + t178 * t3) * s26 - t208) * s34 + (((t26
     & + t125) * s25 + t227) * s25 + t53) * s25 + ddreal(30) * t4 + s26 
     &* (t75 * t84 + t121 + t122) + t119 * t123 * t84) * s34 + s16 * (t1
     &6 * (-t102 + t32) - t100) - t34 * (s25 * t74 + (s26 * t120 + t105 
     &+ t106) * s26) + t57 * t88) * s34 + t67 * s56) * MB1110(1) * t78) 
     &/ ddreal(12) + t6 / ddreal(6)

           intHs16s25s26s34s56x1114D6eps0 = result
       end function intHs16s25s26s34s56x1114D6eps0

       function intHs16s25s26s34s56x1114D6eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1114D6eps1
           type(dd_complex) ::  t1,t2,t3

           type(dd_complex) :: result

      t1 = s16 + s26 + s56
      t2 = s16 + s26 - s34 + s56
      t3 = ddreal(1) / (t1)**2
      result = ddreal(1) / s16 * t3 * (s34 / ddreal(6) - t1 / ddreal(4) 
     &+ (t2)**2 * MB1011(1) / ddreal(6))

           intHs16s25s26s34s56x1114D6eps1 = result
       end function intHs16s25s26s34s56x1114D6eps1

       function intHs16s25s26s34s56x1120D2eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1120D2eps0
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s16 + s25 + s26 - s34 + s56
      t1 = ddreal(1) / t1
      result = ddreal(2) * s25 * MB1110(1) * t1 - (-s25 * MB1110(0) + (s
     &16 + s26 + s56) * MB1011(1) - MB1001(0)) * t1

           intHs16s25s26s34s56x1120D2eps0 = result
       end function intHs16s25s26s34s56x1120D2eps0

       function intHs16s25s26s34s56x1120D2eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1120D2eps1
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s16 + s25 + s26 - s34 + s56
      t1 = ddreal(1) / t1
      result = t1 * (s25 * MB1110(1) + ddreal(1))

           intHs16s25s26s34s56x1120D2eps1 = result
       end function intHs16s25s26s34s56x1120D2eps1

       function intHs16s25s26s34s56x1121D2eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1121D2eps0
           type(dd_complex) ::  t1,t2,t3

           type(dd_complex) :: result

      t1 = s16 + s26 - s34 + s56
      t2 = s16 - s25
      t1 = ddreal(1) / t1
      t3 = ddreal(1) / s16
      result = ddreal(2) * s25 * s34 * MB1101(0) * t3 * (t1)**2 + t1 * (
     &(s25 * MB1110(0) + s34 * MB1011(0) - (-s16 * t2 - s26 * t2 + (ddre
     &al(2) * s16 + s25 + s26 - s34) * s34 - (s16 - s25 - s34) * s56) * 
     &MB1111D2(0) * t1) * t3 - I300s16(0))

           intHs16s25s26s34s56x1121D2eps0 = result
       end function intHs16s25s26s34s56x1121D2eps0

       function intHs16s25s26s34s56x1121D2eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1121D2eps1
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s16 + s26 - s34 + s56
      t1 = ddreal(1) / t1
      result = t1 * ((s25 * MB1110(1) + s34 * MB1011(1)) / s16 - I300s16
     &(1))

           intHs16s25s26s34s56x1121D2eps1 = result
       end function intHs16s25s26s34s56x1121D2eps1

       function intHs16s25s26s34s56x1121D4eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1121D4eps0
           type(dd_complex) ::  t1,t2,t3,t4,t5

           type(dd_complex) :: result

      t1 = s16 + s26 - s34 + s56
      t2 = ddreal(2)
      t3 = s16 * s25 + (-s34 + s16 + s25 + s26) * s26 + (s26 * t2 + s16 
     &+ s25 - s34 + s56) * s56
      t4 = -s16 + s25
      t5 = ddreal(1) / t3
      t3 = ddreal(1) / t3
      result = -t1 * t3 - t5 * (t1 * (s16 * I300s16(1) - s25 * MB1110(1)
     & + t1 * MB1011(1) - MB1001(0)) - (s16 * t4 + s26 * t4 + (s16 * t2 
     &+ s25 + s26 - s34) * s34 + (-s16 + s25 + s34) * s56) * MB1111D2(0)
     & + s25 * s34 * MB1101(0)) / ddreal(2)

           intHs16s25s26s34s56x1121D4eps0 = result
       end function intHs16s25s26s34s56x1121D4eps0

       function intHs16s25s26s34s56x1121D4eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1121D4eps1

           type(dd_complex) :: result

      result = ddreal(0)

           intHs16s25s26s34s56x1121D4eps1 = result
       end function intHs16s25s26s34s56x1121D4eps1

       function intHs16s25s26s34s56x1122D4eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1122D4eps0
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t2,t3,t4,t5
           type(dd_complex) ::  t6,t7,t8,t9

           type(dd_complex) :: result

      t1 = s16 + s25 + s26 - s34 + s56
      t2 = ddreal(2)
      t3 = t2 * s26
      t4 = s16 * s25 + (-s34 + s16 + s25 + s26) * s26 + (t3 + s16 + s25 
     &- s34 + s56) * s56
      t5 = s16 + s26 - s34 + s56
      t6 = s16 + s26 + s56
      t7 = s25 + s26
      t8 = t2 * s16
      t9 = s16 + s26
      t10 = s16 * t9
      t11 = (s26)**2
      t12 = (s16)**2
      t13 = (t11 + t12) * s16
      t14 = ddreal(3) * s16
      t6 = ddreal(1) / t6
      t15 = ddreal(1) / t5
      t16 = ddreal(1) / s16
      t4 = ddreal(1) / t4
      result = t2 * t1 * t16 * t4 - t6 * (t4 * (t5 * MB1001(0) - (t3 * t
     &12 + ((s16 - s26) * s34 - (t8 - t7) * t9) * s34 + ((s16 + s34) * s
     &56 + t10 * t2 - s34 * (s34 - t3 + s16 - s25)) * s56 + t13) * MB101
     &1(1) * t16) - s34 * MB1011(0) * t16) + (t2 * (s26 * t12 + ((-t8 - 
     &t7 + s34) * s34 + t10) * s56) + ((t2 * t7 - s34 + t14) * s34 - t2 
     &* s25 * t9 - (ddreal(4) * s26 + t14) * s16 - (s25)**2 - t11) * s34
     & + (s16 - s34) * (s56)**2 + t13) * MB1111D2(0) * t16 * t15 * t4 - 
     &t4 * (s25 * t16 * (-s34 * MB1101(0) * t15 + MB1110(1)) - I300s16(1
     &)) * t1

           intHs16s25s26s34s56x1122D4eps0 = result
       end function intHs16s25s26s34s56x1122D4eps0

       function intHs16s25s26s34s56x1122D4eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1122D4eps1
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s16 + s26 + s56
      t1 = ddreal(1) / t1
      result = ddreal(1) / s16 * t1 * (s34 * MB1011(1) + ddreal(1))

           intHs16s25s26s34s56x1122D4eps1 = result
       end function intHs16s25s26s34s56x1122D4eps1

       function intHs16s25s26s34s56x1123D6eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1123D6eps0
           type(dd_complex) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           type(dd_complex) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           type(dd_complex) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           type(dd_complex) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           type(dd_complex) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           type(dd_complex) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           type(dd_complex) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           type(dd_complex) ::  t175,t176,t177,t178,t179,t18,t19,t2,t20,t21,t22,t23
           type(dd_complex) ::  t24,t25,t26,t27,t28,t29,t3,t30,t31,t32,t33,t34
           type(dd_complex) ::  t35,t36,t37,t38,t39,t4,t40,t41,t42,t43,t44,t45
           type(dd_complex) ::  t46,t47,t48,t49,t5,t50,t51,t52,t53,t54,t55,t56
           type(dd_complex) ::  t57,t58,t59,t6,t60,t61,t62,t63,t64,t65,t66,t67
           type(dd_complex) ::  t68,t69,t7,t70,t71,t72,t73,t74,t75,t76,t77,t78
           type(dd_complex) ::  t79,t8,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89
           type(dd_complex) ::  t9,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99

           type(dd_complex) :: result

      t1 = s16 - s25
      t2 = ddreal(2)
      t3 = t2 * s16
      t4 = s25 + t3
      t5 = s16 + s25
      t6 = s25 * t5
      t7 = ddreal(3)
      t8 = (s16)**2
      t9 = (t8)**2
      t10 = s16 * t8
      t11 = t10 * t9
      t12 = t8 * t9
      t13 = s16 * t9
      t14 = t7 * t8
      t15 = ddreal(4)
      t16 = t15 * s25
      t17 = ddreal(9) * t8
      t18 = (s25)**2
      t19 = (t18)**2
      t20 = s25 * t18
      t21 = t8 + t18
      t22 = s25 * t21
      t23 = ddreal(9) * t22
      t24 = (ddreal(11) * t18 + ddreal(24) * t8) * s16
      t25 = ddreal(11) * s16
      t26 = ddreal(17) * s25
      t27 = (t25 + t26) * s25
      t28 = ddreal(27) * t8
      t29 = t27 + t28
      t30 = s16 * s25
      t31 = t30 * t21
      t32 = ddreal(25)
      t33 = ddreal(7)
      t34 = ddreal(5)
      t35 = t15 * s26
      t36 = t33 * t8
      t37 = t18 * t34
      t38 = ddreal(66) * s16
      t39 = ddreal(42) * s25
      t40 = (t38 + t39) * s25
      t41 = ddreal(72) * t8
      t42 = ddreal(27) * s16
      t43 = ddreal(22) * s25
      t44 = (s26)**2
      t45 = (t44)**2
      t46 = s26 * t44
      t47 = ddreal(20) * s25
      t48 = t2 * t46
      t49 = ddreal(48) * s16
      t50 = ddreal(9) * s26
      t51 = ddreal(57) * s16
      t52 = ddreal(50) * t8
      t53 = ddreal(12) * s26
      t54 = t32 * s16
      t55 = ddreal(6)
      t56 = t55 * s26
      t57 = ddreal(44)
      t58 = ddreal(24) * s16
      t59 = (s34)**2
      t60 = (t59)**2
      t61 = s34 * t60
      t62 = s34 * t59
      t63 = t1 * ((-t3 - t16) * s25 + t17)
      t64 = t9 + t60
      t65 = ddreal(10)
      t66 = (t1)**2
      t67 = ddreal(8) * s16
      t68 = t65 * t8
      t69 = s26 * t5
      t70 = t10 - t62
      t71 = ddreal(11) * t8
      t72 = s34 * t5
      t73 = -t2 * (t30 + t72) + t8 + t59 + t18
      t74 = (s56)**2
      t75 = t1 * t44
      t76 = t1 * s26
      t77 = t7 * s16
      t78 = (t77 + t16) * s25
      t23 = t2 * (t73 * s56 * t74 + s16 * t19 + t76 * (t4 * (-t6 + t14) 
     &+ t75)) + t34 * (t13 - t61) + ((((t47 + t53 + t54) * s34 - (ddreal
     &(30) * s25 + t51) * s25 - (ddreal(38) * s25 + t49 + t50) * s26 - t
     &52) * s34 + ((ddreal(41) * s16 + t47) * s25 + ddreal(51) * t8) * s
     &25 + ddreal(50) * t10 + s26 * ((t42 + t43) * s26 + t40 + t41) + t4
     &8) * s34 - ((t35 * t5 + t29) * s26 + t2 * (t23 + t24)) * s26 - t18
     & * (t36 + t37) - t32 * t9 - ddreal(11) * t31) * s34 + ((s25 * ((-t
     &3 + t16) * s25 - t71) - s34 * (-(t42 + t43 + t56) * s34 + t27 + t2
     &8 + ddreal(12) * t69) + ddreal(9) * t70 + t56 * t66) * s56 - t2 * 
     &(((ddreal(19) * s25 + t50 + t58) * t59 + (t56 * t5 + t29) * s26 + 
     &t23 + t24) * s34 - t19 - t63 * s26) + ((s25 * t57 + ddreal(54) * s
     &16 + t56) * s26 + t40 + t41) * t59 + ddreal(12) * t64 + t30 * ((-t
     &67 + t16) * s25 - t68) + t55 * t66 * t44) * s56 + t63 * t44 - t78 
     &* t10
      t24 = t2 * s26
      t27 = (-s34 + t5 + s26) * s26 + (s16 + s25 - s34 + t24 + s56) * s5
     &6 + t30
      t29 = t7 * s25
      t40 = t34 * s16
      t41 = -t18 + t14
      t63 = -s25 + t3
      t79 = s16 + s26
      t80 = ddreal(11) * s25
      t81 = ddreal(13) * t8
      t82 = t34 * s25
      t83 = t15 * s16
      t84 = t82 + t83
      t85 = ddreal(19) * s16
      t86 = ddreal(18) * s25
      t87 = t33 * s16
      t88 = ddreal(15) * t8
      t89 = t65 * t10
      t90 = ddreal(32) * s16
      t91 = ddreal(9) * s25
      t92 = (t90 + t91) * s25
      t93 = ddreal(38) * t8
      t94 = s16 * t18
      t95 = t46 + t94
      t96 = ddreal(15) * s25
      t97 = t65 * s16
      t98 = ddreal(8) * s26
      t99 = ddreal(17) * s16
      t100 = t7 * s26
      t101 = s16 + t100
      t102 = t18 * t5
      t103 = t8 * (t87 - t91) + t102
      t104 = ddreal(33) * s25
      t105 = t65 * s25
      t106 = ddreal(8) * t8
      t107 = t55 * s25
      t108 = ddreal(12) * s16
      t109 = ddreal(17) * t8
      t110 = ddreal(34) * t8
      t111 = ddreal(42) * t10
      t112 = ddreal(16)
      t113 = t112 * s26
      t114 = -t8 + t18
      t115 = t8 * s25
      t116 = t10 * t114
      t117 = t15 * t44
      t118 = t55 * s16
      t119 = ddreal(20) * s16
      t120 = s25 * s26
      t121 = ddreal(21) * s16
      t122 = s25 + s26
      t123 = t94 * t1
      t124 = t10 * t1
      t125 = t35 * t1
      t126 = t8 + t59
      t127 = t2 * t126
      t128 = t83 * s34
      t129 = -t10 - t20
      t130 = (t79)**2
      t131 = (-((-s34 * t34 + s25 - t77) * s25 + t127 - t128) * s56 - t3
     &3 * t70 + (-(t98 + t26 + t121) * s34 + s16 * (t113 + t121) + (ddre
     &al(20) * s26 + t67 + t80) * s25) * s34 - t102 - t125 * t63 + ddrea
     &l(9) * t115) * s56 + (-ddreal(8) * t62 + ((t58 + t104) * s25 + ddr
     &eal(63) * t8) * s26 + ddreal(35) * t10 + t18 * (t119 + t82) + t55 
     &* (t44 * t84 + t115)) * s34 - t7 * (s26 * t103 + t59 * (t121 * s26
     & - (t122 * t33 + t25) * s34 + t18 * t55 + t117 + ddreal(17) * t120
     & + ddreal(12) * t30 + ddreal(17) * t8) - t123) - ddreal(9) * t124 
     &- t75 * t55 * t63
      t38 = t131 * s56 + t34 * t116 - t7 * (t103 * t44 + t115 * t114 - t
     &61) - ((((t113 + t91 + t99) * s34 - (ddreal(21) * s26 + t38 + t39)
     & * s26 - t92 - t93) * s34 + ((t29 + t54) * s25 + t110) * s25 + ((d
     &dreal(63) * s16 + ddreal(51) * s25 + t98) * s26 + t55 * ((t107 + t
     &108) * s25 + t109)) * s26 + t111) * s34 - t79 * (((t87 + t105) * s
     &25 + t106) * s25 + ((t83 + t104) * s25 + ddreal(47) * t8 + t35 * t
     &84) * s26 + ddreal(23) * t10)) * s34 - t76 * (t117 * t63 + t118 * 
     &t41)
      t39 = t129 * t10
      t103 = t76 * t8
      t131 = t2 * t18
      t132 = s16 + s26 + s56
      t133 = t34 * t8
      t134 = t7 * t18
      t135 = t2 * s25
      t136 = (-s16 - t135) * s25
      t137 = t15 * t8 + t136
      t138 = ddreal(20) * t8
      t139 = t18 * t65
      t140 = ddreal(8) * s25
      t141 = (t40 + t140) * s25
      t142 = ddreal(12) * t8
      t143 = t142 + t141
      t144 = (ddreal(28) * s16 + t86) * s25
      t145 = ddreal(30) * t8
      t146 = ddreal(21) * t8
      t147 = ddreal(20) * t10
      t148 = t2 * t44
      t149 = s16 * t79
      t150 = ddreal(23) * s16
      t151 = ddreal(12) * s25
      t152 = t34 * s26
      t153 = t7 * t5
      t154 = s16 * t20
      t155 = t7 * t44
      t70 = t15 * t70
      t156 = t76 * ((t133 - t134) * s16 - s25 * t114)
      t19 = -t2 * (t61 + t116) + ((((t152 + t140 + t97) * s34 - (t150 + 
     &t113 + t151) * s25 - t117 - ddreal(20) * t149) * s34 + ((t140 + t9
     &9) * s25 + t146) * s25 + s26 * (t145 + t144 + t44) + t147 + t148 *
     & (t82 + t118)) * s34 - t2 * (t46 * t5 + t19) - t31 * t34 - t8 * (t
     &68 + t134) - s26 * (t143 * s26 + (t138 + t139) * s16 + ddreal(8) *
     & t22)) * s34 + ((s56 * t73 + (-t133 + t131) * s25 - (-(t100 + t105
     & + t108) * s34 + t55 * t69 + t141 + t142) * s34 - t94 + t70 + t100
     & * t66) * s56 - t15 * (t62 * (t24 + t40 + t16) + t115 * t5) - t2 *
     & (((t68 + t37) * s16 + (t153 * s26 + t143) * s26 + t15 * t22) * s3
     &4 - t154 - t76 * t137) + t34 * t64 + t59 * ((t58 + t47 + t100) * s
     &26 + t144 + t145) + t19 + t155 * t66) * s56 + t46 * t66 + t30 * (-
     &t10 + t20) + t75 * t137 + t156
      t31 = (s16 - t107) * s25
      t37 = t133 + t136
      t47 = (-t135 - t40) * s25 + t71
      t64 = -t135 + t77
      t66 = t112 * s25
      t68 = t33 * s25
      t105 = t118 + t68
      t113 = ddreal(31) * s16
      t116 = t32 * s25
      t117 = s25 + t83
      t136 = ddreal(29)
      t137 = t136 * s25
      t141 = ddreal(55) * s16
      t143 = (t141 + t116) * s25
      t144 = ddreal(78) * t8
      t156 = ddreal(93) * s16
      t157 = ddreal(72) * s25
      t158 = ddreal(37) * s16
      t159 = ddreal(15) * s26
      t160 = t64 * t1
      t161 = t160 * t46
      t162 = ddreal(18) * s26
      t163 = ddreal(30) * s26
      t164 = t8 * t18
      t165 = t75 * t64
      t166 = ddreal(15) * s16
      t64 = t76 * t64
      t126 = (((s34 * t33 - t135 + t40) * s25 - t126 * t7 + t118 * s34) 
     &* s56 + ((-t135 - t77) * s25 + t112 * t8) * s25 + t15 * ((s16 * (t
     &67 + t56) + (s26 * t33 + t16 + t77) * s25) * s34 - t64) + t59 * (s
     &34 * t65 - ddreal(24) * s25 - t113 - t53) - ddreal(11) * t10) * s5
     &6 + t10 * (-t166 + t86) - t55 * (t154 + t165) - t7 * (t76 * t47 - 
     &t164) - ddreal(11) * t60 + (((t49 + t163 + t137) * s34 - (t162 + t
     &156 + t157) * s26 - t143 - t144) * s34 + ((t113 + t68) * s25 + t10
     &6) * s25 + (t56 * t105 + ddreal(12) * t106 + ddreal(12) * t78) * s
     &26 + ddreal(56) * t10) * s34
      t53 = t126 * s56 + t10 * ((t67 + t68) * s25 - t17) - t15 * (-t61 +
     & t161) + t2 * s34 * (t59 * ((t58 + t107) * s25 + (t49 + t159 + t13
     &7) * s26 + t145) + t79 * (t102 * t33 + t15 * (((s16 + t107) * s25 
     &+ t17) * s26 + t115) + t147 + t148 * t105)) - t55 * s16 * (t76 * t
     &37 + t154) - t59 * (t59 * (ddreal(22) * s26 + t151 + t54) + ((t16 
     &+ t158) * s25 + ddreal(52) * t8) * s25 + ((t156 + t157 + t53) * s2
     &6 + t2 * (t144 + t143)) * s26 + ddreal(70) * t10) - t75 * t7 * t47
      t58 = t160 * t45
      t67 = t10 * t18
      t17 = t2 * (t61 * (s16 + t24) + t39) + t7 * s16 * (-t75 * t37 + t6
     &7) + s25 * t13 + s34 * (((-(s16 * (t140 + t97) + (ddreal(11) * s26
     & + t151 + t54) * s26) * s34 + (ddreal(12) * s25 * t117 + (t49 + t1
     &37) * s26 + ddreal(12) * t133) * s26 + t65 * t95 + t8 * (ddreal(23
     &) * s25 + t119)) * s34 - t79 * (s16 * ((t151 + t121) * s25 + t138)
     & + (t155 + (t113 + t116) * s25 + t52) * s26 + t15 * (t44 * (t107 +
     & t87) + t20))) * s34 + t130 * ((t18 * t33 + t133) * s25 + t89 + (t
     &105 * s26 + (-t3 + t66) * s25 + t138) * s26 - t94)) + t53 * s56 - 
     &t58 - t1 * t47 * t46 - t103 * (t31 + t17)
      t37 = -t29 + t83
      t47 = t135 + t77
      t49 = ddreal(26)
      t52 = ddreal(8) * t18
      t53 = t49 * t10
      t54 = t8 * ((t52 - t110) * s25 + t53)
      t68 = ddreal(14)
      t102 = t68 * t8
      t105 = (-t135 - t87) * s25 + t102
      t113 = ddreal(24) * t10
      t14 = (s16 * t68 - t135) * s25 + t14
      t116 = t14 * s26
      t119 = t34 * t10
      t82 = -t77 - t82
      t121 = t2 * t8
      t25 = t107 + t25
      t126 = -t108 - t137
      t138 = s16 * ((-t29 + t141) * s25 + ddreal(105) * t8)
      t141 = t15 * t10
      t143 = ddreal(62)
      t144 = (s16 * t57 - t107) * s25
      t145 = t152 * t14
      t147 = -ddreal(45) * s25
      t152 = ddreal(75) * t8
      t156 = -t82
      t157 = (t87 - t104) * s25
      t57 = t57 * t8
      t167 = ddreal(88) * s25
      t168 = ddreal(52) * s16
      t169 = t167 + t168
      t170 = (t57 + t157) * t7
      t171 = ddreal(36)
      t172 = t171 * s26
      t173 = t7 * (s16 * (t29 + t108) + (t3 - t140 - t50) * s26) * t61
      t174 = ddreal(13) * t10
      t175 = (t42 - t26) * s25 + t88
      t176 = t156 * t65 * s26
      t177 = s16 * s26
      t178 = ((s16 * t49 - t16) * s25 + ddreal(28) * t8) * s25
      t179 = ddreal(33) * t10
      t14 = ((t7 * (t10 + t62) + s34 * (s34 * t82 - t14) + t30 * (t135 -
     & t40)) * s56 + s16 * ((-t146 + t131) * s25 + t10 * t68) + t34 * s1
     &6 * (t64 + t94) - ddreal(9) * t60 - ((-(ddreal(13) * s16 + t43 + t
     &159) * s34 - (-s26 * t32 - t26 + t42) * s25 + ddreal(15) * s16 * (
     &-s16 + s26)) * s34 + t145 + t178 + t179) * s34) * s56 + s16 * (s16
     & * ((t52 - t110) * s25 + t53) + t165 * t65 + t125 * t105) + ddreal
     &(9) * t61 + (((-(t108 + t137 + t172) * s34 + (t167 + t168 + t163) 
     &* s26 - t157 - t57) * s34 + ((t150 - t96) * s25 + t152) * s25 + (t
     &15 * t175 - t176) * s26 + ddreal(114) * t10) * s34 + ((t131 - ddre
     &al(82) * t8) * s25 - t119) * s25 - t15 * ((t179 + t178) * s26 + t1
     &54) - ddreal(93) * t9 - t14 * t65 * t44) * s34
      t14 = t14 * s56 + t10 * ((-t49 * t8 + ddreal(12) * t18) * s25 + t1
     &13) + t55 * t177 * (s16 * ((t15 * t18 - t109) * s25 + t174) + t76 
     &* t105) - t65 * s16 * (-t161 + t67) + t7 * ((-s16 + t16 + t50 - s3
     &4) * s34 + (-t162 - t108 - t137) * s26 + t31 + ddreal(22) * t8) * 
     &t60 + (((((ddreal(9) * s16 + t151) * s25 - ddreal(81) * t8) * s25 
     &+ ((ddreal(3) / ddreal(2) * t169 + t163) * s26 - t170) * s26 - ddr
     &eal(174) * t10) * s34 + (((-t166 - t29) * s25 + ddreal(99) * t8) *
     & s25 + ddreal(79) * t10) * s25 + ((t175 * t55 - t176) * s26 + ((dd
     &real(69) * s16 + t147) * s25 + ddreal(225) * t8) * s25 + ddreal(34
     &2) * t10) * s26 + ddreal(201) * t9) * s34 - t79 * ((((s16 * t171 -
     & t107) * s25 + ddreal(70) * t8) * s25 - t174) * s25 + (t15 * ((t14
     &4 + t36) * s25 + t111) + t116 * t65) * s26 + ddreal(111) * t9)) * 
     &s34
      t13 = t14 * s56 + t15 * t76 * s16 * (t8 * t37 * t47 + t105 * t44) 
     &+ ddreal(11) * t11 + s34 * ((((t59 * (t40 - t56) + ((t126 * t7 - t
     &172) * s26 + t63 * t25 * t55) * s26 + t138) * s34 + s16 * (((-ddre
     &al(40) * s16 + ddreal(13) * s25) * s25 - ddreal(102) * t8) * s25 -
     & ddreal(160) * t10) + (((t169 + t159) * s26 - t170) * s26 - t55 * 
     &((-t78 + t28) * s25 + ddreal(58) * t10)) * s26) * s34 + t79 * ((((
     &t166 - t107) * s25 + ddreal(61) * t8) * s25 + ddreal(66) * t10) * 
     &s25 + ddreal(135) * t9 + s26 * (((ddreal(133) * s16 - ddreal(68) *
     & s25) * s25 + t152) * s26 + ((ddreal(137) * s16 + t147) * s25 + dd
     &real(92) * t8) * s25 + ddreal(267) * t10) - t34 * t156 * t46)) * s
     &34 - t130 * (((t144 + t106) * s25 + t10) * s25 + (t2 * (((s16 * t1
     &43 - t140) * s25 - t102) * s25 + ddreal(51) * t10) + t145) * s26 +
     & ddreal(60) * t9)) + t55 * t54 * t44 / ddreal(2) + ddreal(8) * t9 
     &* t20 - t139 * t13 - t173 - t91 * t12 + t40 * t58
      t14 = t76 * t9 * ((t3 - t140) * s25 + t71)
      t11 = -s25 * t11 + t13 * s56 + t46 * t54 + t60 * (t10 * (ddreal(30
     &) * s16 + ddreal(35) * s25) + s26 * (t126 * t44 + t138) - ddreal(9
     &) * t45 + ddreal(9) * t164 + t155 * t63 * t25) + t160 * s16 * s26 
     &* t45 + s16 * t1 * t105 * t45 + t130 * (t18 * ((t166 - t29) * s25 
     &+ t121) + ((t82 * s26 + (t158 - t26) * s25 + t146) * s26 + ((t51 -
     & t96) * s25 + t55 * t8) * s25 + ddreal(75) * t10) * s26 + ddreal(3
     &0) * t10 * t5) * t59 - t79 * t130 * ((((s16 * t112 - t135) * s25 -
     & t106) * s25 + t119) * s25 + ddreal(12) * t9 + s26 * (((-t16 + t90
     &) * s25 - t102) * s25 + t113 + t116)) * s34 + (-s26 + t3) * t101 *
     & t59 * t60 - t79 * ((t15 * (((-t29 + t118) * s25 + t81) * s25 + dd
     &real(30) * t10) + t44 * (-t2 * (t40 + t80) - t100)) * s26 + t8 * (
     &(ddreal(50) * s16 + t66) * s25 + ddreal(40) * t8) + t44 * ((s16 * 
     &t136 - t104) * s25 + ddreal(54) * t8) - t154) * t62 + t14
      t13 = t1 * s25
      t14 = t55 * t21
      t25 = t30 * t34
      t28 = t2 * t5
      t31 = (-t8 + t13) * s26
      t32 = s16 + s25 + s26 - s34 + s56
      t42 = ddreal(1) / s16
      t43 = ddreal(1) / t73
      t27 = ddreal(1) / (t27)**2
      t49 = ddreal(1) / (t132)**2
      t1 = t17 * MB1001(0) + (-t75 * t41 * t7 + t131 * t10) * s16 + s34 
     &* ((((t101 * s34 - t34 * s16 * t5 - (t98 + t91 + t99) * s26) * s34
     & + t33 * t95 + t8 * (t96 + t97) + s26 * ((ddreal(33) * s16 + ddrea
     &l(21) * s25) * s26 + t92 + t93)) * s34 - t79 * (((t29 + t87) * s25
     & + t88) * s25 + t48 + s26 * ((t85 + t26) * s26 + (t85 + t86) * s25
     & + ddreal(32) * t8) + t89)) * s34 + t130 * (t34 * (t10 + t22) + s2
     &6 * (t84 * s26 + (-t3 + t80) * s25 + t81) - t3 * t18)) + t38 * s56
     & + t39 - t1 * t63 * t45 - t1 * (-s25 * t4 + t36) * t46 - t103 * (-
     &t29 + t40) * t5
      t1 = t43 * t27 * (t19 * I300s16(1) + (-t32 * s25 * (-t2 * (t62 * (
     &t24 + t29 + t83 - s34) - t124) + (((s25 - t77) * s25 + t121) * s26
     & + (t18 - t133) * s25 + t141) * s26 + (((t107 + t97 + t50) * s25 +
     & t148 + ddreal(12) * t149) * s34 - (t106 + t134) * s16 - t2 * t22 
     &- s26 * ((t29 + t83) * s26 + (t83 + t107) * s25 + t142)) * s34 + (
     &(-t7 * s25 * (s16 + s34) + t127 - t128 + t18) * s56 - t2 * ((t7 * 
     &s25 * t122 + (t135 + t35 + t118) * s16) * s34 - t76 * t63) + t70 +
     & t59 * (t35 + t91 + t108) + t20 - t115 * t34) * s56 - t123) * MB11
     &10(1) + t23 + (t2 * (t75 * t37 * t47 - t129 * t8) * t10 - t7 * (t1
     &2 * t18 + t61 * (t44 * (s16 - t16) + t7 * ((s16 * t117 - t44) * s2
     &6 + t115) + t141)) + t11) * MB1011(1) * t49) * t42 + t49 * t1)
      result = t1 / ddreal(4) + t42 * (t27 * (s25 * s34 * (-(-t5 * (t2 *
     & t21 - t30 * t7) + t31) * s26 - ((-t55 * (t8 + t18 + t120 + t177) 
     &- t44 - ddreal(8) * t30 - t59) * s34 + (t28 * s26 + t14 + t25) * s
     &26 + t15 * ((t8 + t6) * s25 + t10)) * s34 - (-(-t2 * t72 - t13 + t
     &59 + t8) * s56 + t2 * (t59 * (-t153 - s26 + s34) - t10 - t20 + t31
     &) + (t15 * t69 + t14 + t25) * s34 + t30 * t5) * s56 + t114 * t18 +
     & t9 - t2 * (s26 + t28) * t62) * MB1101(0) * t43 + t32 * (-t2 * (-s
     &26 * t8 + ((t122 + t3 - s34) * s34 - t149) * s56) - (-t8 - t44) * 
     &s16 - ((-t122 * t2 + s34 - t77) * s34 + t2 * s25 * t79 + (t35 + t7
     &7) * s16 + t18 + t44) * s34 + t74 * (s16 - s34)) * MB1111D2(0)) + 
     &s34 * (s16 + s26 - s34 + s56) * MB1011(0) * t49) / ddreal(2)

           intHs16s25s26s34s56x1123D6eps0 = result
       end function intHs16s25s26s34s56x1123D6eps0

       function intHs16s25s26s34s56x1123D6eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1123D6eps1
           type(dd_complex) ::  t1,t2

           type(dd_complex) :: result

      t1 = s16 + s26 + s56
      t2 = ddreal(1) / (t1)**2
      result = ddreal(1) / s16 * t2 * (-s34 / ddreal(2) + t1 / ddreal(4)
     & + s34 * (s16 + s26 - s34 + s56) * MB1011(1) / ddreal(2))

           intHs16s25s26s34s56x1123D6eps1 = result
       end function intHs16s25s26s34s56x1123D6eps1

       function intHs16s25s26s34s56x1130D4eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1130D4eps0
           type(dd_complex) ::  t1,t2

           type(dd_complex) :: result

      t1 = s16 + s25 + s26 - s34 + s56
      t2 = ddreal(3)
      t1 = ddreal(1) / (t1)**2
      result = -(s25 * t2 + s16 + s26 - s34 + s56) * t1 * ((s16 + s26 + 
     &s56) * MB1011(1) - MB1001(0)) / ddreal(4) + s25 * t1 * (s25 * MB11
     &10(0) + ddreal(1)) / ddreal(2) + ddreal(3) / ddreal(2) * (s25)**2 
     &* MB1110(1) * t1

           intHs16s25s26s34s56x1130D4eps0 = result
       end function intHs16s25s26s34s56x1130D4eps0

       function intHs16s25s26s34s56x1130D4eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1130D4eps1
           type(dd_complex) ::  t1,t2,t3

           type(dd_complex) :: result

      t1 = ddreal(3)
      t2 = s16 + s25 + s26 - s34 + s56
      t3 = ddreal(1) / ddreal(4)
      t2 = ddreal(1) / (t2)**2
      result = t2 * (t3 * (s25 * t1 + s16 + s26 - s34 + s56) + (s25)**2 
     &* MB1110(1) / ddreal(2))

           intHs16s25s26s34s56x1130D4eps1 = result
       end function intHs16s25s26s34s56x1130D4eps1

       function intHs16s25s26s34s56x1131D4eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1131D4eps0
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           type(dd_complex) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           type(dd_complex) ::  t4,t5,t6,t7,t8,t9

           type(dd_complex) :: result

      t1 = s16 + s25 - s34
      t2 = ddreal(2)
      t3 = -t2 * (s26 + s56) - t1
      t4 = s16 + s25
      t5 = t2 * s26
      t6 = s16 * s25 + (-s34 + t4 + s26) * s26 + (t5 + t1 + s56) * s56
      t7 = s16 + s26 - s34 + s56
      t8 = s16 + s26
      t9 = (s16)**2
      t10 = (s26)**2
      t11 = (s56)**2
      t12 = s16 + s26 + s56
      t13 = s16 + s25 + s26 - s34 + s56
      t14 = s16 * t8
      t15 = -s16 + s26 + s25
      t16 = ddreal(3)
      t17 = t16 * t14
      t18 = s16 - s25
      t14 = ddreal(6) * t14
      t19 = (s34)**2
      t20 = t16 * s16
      t21 = s25 + s26
      t22 = t2 * s16
      t23 = t22 + t21
      t24 = (s25)**2
      t25 = t24 + t9
      t26 = (t18)**2
      t27 = t26 * s26
      t28 = t9 * s25
      t7 = ddreal(1) / t7
      t29 = ddreal(1) / t6
      t30 = ddreal(1) / s16
      t6 = ddreal(1) / t6
      t13 = ddreal(1) / t13
      t12 = ddreal(1) / t12
      t21 = -(t2 * (s34 * t19 * t23 + (s34 * t8 * (t22 - s25 + s26) + t2
     &8 - t27) * s16 - (-t28 * t2 + s16 * t25 - (s16 * (-s25 + t20 + t5)
     & + (s34 - t20 - t21) * s34) * s34 + t27) * s56) - t10 * t26 - t11 
     &* (-t2 * s16 * (s25 + s34) + t19 + t24 + t9) - t19 * (t2 * s25 * t
     &8 + t10 + t14 + t19 + t24) - t25 * t9) * MB1111D2(0) * t30 * t7 * 
     &t6 + I300s16(0)
      t22 = t19 * MB1011(0) * t30 * t12 * t7
      t1 = t13 * (-t24 * MB1110(0) * t7 * t30 - s25 * (s26 * t4 + (-t2 *
     & t4 - s26 + s34) * s34 + s56 * t1 + t24 + t9) * MB1110(1) * t30 * 
     &t6 + t6 * (-(t2 * (s16 * s26 + (s16 + s26 - s34) * s56) + (-t2 * t
     &8 - s25 + s34) * s34 + t10 + t11 + t9) * MB1001(0) + (t16 * t9 * s
     &26 * t8 + s16 * (s16 * t9 + s26 * t10) - ((-t2 * s26 * t18 + (s16 
     &- s26) * s34 - t16 * t9 + t10) * s34 + t8 * (-t15 * s25 + t17)) * 
     &s34 - ((-s16 * s56 - t17 + s34 * (s34 - s25 + t20)) * s56 - t16 * 
     &s16 * (t10 + t9) - s34 * (t19 + (s25 + t5) * s25 - t14) - ddreal(6
     &) * t9 * s26 + t2 * t15 * t19) * s56) * MB1011(1) * t30) * t12) + 
     &t7 * t21 + s25 * s34 * (s16 * t18 + s26 * t18 + (-t23 + s34) * s34
     & + (s16 - s25 - s34) * s56) * MB1101(0) * t30 * (t7)**2 * t6 - t3 
     &* I300s16(1) * t29 - t22
      result = -t1 / ddreal(2) + t3 * t30 * t6

           intHs16s25s26s34s56x1131D4eps0 = result
       end function intHs16s25s26s34s56x1131D4eps0

       function intHs16s25s26s34s56x1131D4eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1131D4eps1
           type(dd_complex) ::  t1,t2,t3,t4,t5,t6,t7

           type(dd_complex) :: result

      t1 = s16 + s26 + s56
      t2 = ddreal(2)
      t3 = s16 + s25 + s26 - s34 + s56
      t4 = s16 + s26 - s34 + s56
      t5 = ddreal(1) / t1
      t3 = ddreal(1) / t3
      t6 = ddreal(1) / s16
      t4 = ddreal(1) / t4
      t7 = -ddreal(1) / ddreal(2)
      result = t7 * (t4 * (-((s34)**2 * MB1011(1) * t5 + (s25)**2 * MB11
     &10(1) * t3) * t6 + I300s16(1)) + (t1 * t2 + s25 - s34) * t6 * t5 *
     & t3)

           intHs16s25s26s34s56x1131D4eps1 = result
       end function intHs16s25s26s34s56x1131D4eps1

       function intHs16s25s26s34s56x1132D6eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1132D6eps0
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           type(dd_complex) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           type(dd_complex) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           type(dd_complex) ::  t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52
           type(dd_complex) ::  t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63
           type(dd_complex) ::  t64,t65,t66,t67,t68,t69,t7,t70,t71,t72,t73,t74
           type(dd_complex) ::  t8,t9

           type(dd_complex) :: result

      t1 = s16 + s26 - s34 + s56
      t2 = s16 + s26
      t3 = ddreal(4)
      t4 = ddreal(2)
      t5 = t4 * s16
      t6 = t3 * s26
      t7 = ddreal(3) * s25
      t8 = ddreal(3) * s26
      t9 = s16 + t8
      t10 = ddreal(8)
      t11 = ddreal(6) * s16
      t12 = t10 * s26
      t13 = (s16)**2
      t14 = (t13)**2
      t15 = s16 * t13
      t16 = (s26)**2
      t17 = s26 * t16
      t18 = (s34)**2
      t19 = s34 * t18
      t20 = s16 * s26
      t21 = s16 + s26 + s56
      t22 = s16 + s25
      t23 = t4 * s26
      t24 = s16 * s25
      t25 = (-s34 + t22 + s26) * s26 + (t23 + s16 - s34 + s25 + s56) * s
     &56 + t24
      t26 = ddreal(7)
      t27 = t4 * s25
      t28 = t10 * s16
      t29 = t26 * t13
      t30 = ddreal(5)
      t31 = t30 * s25
      t32 = s26 + s25
      t33 = s25 * t32
      t34 = ddreal(14)
      t35 = ddreal(15) * s16
      t36 = (s25)**2
      t37 = (s16 + t31) * s26
      t38 = s16 - s34
      t39 = ddreal(9) * s26
      t40 = -t19 + t15
      t41 = (s56)**2
      t42 = (t41)**2
      t43 = s56 * t41
      t44 = t13 * s25
      t45 = s16 - s25
      t46 = s16 * t45
      t47 = s16 * t2
      t48 = ddreal(6) * t47
      t49 = ddreal(3) * s16
      t50 = t3 * s25
      t51 = s34 * t4
      t52 = t30 * s16
      t53 = t3 * s16
      t54 = t30 * s26
      t55 = t4 * s25 * t2
      t56 = t28 - t7
      t57 = s16 * t26 - t27
      t58 = t36 + t16
      t59 = ddreal(6) * s26
      t60 = t27 * s26
      t61 = s26 * t32
      t62 = -s25 + t11
      t63 = ddreal(32) * t2
      t64 = t7 * s26
      t65 = ddreal(24)
      t66 = (t2)**2
      t67 = s16 * t66
      t68 = s25 * t22
      t69 = t62 * s26
      t70 = ddreal(6) * s25
      t71 = ddreal(3) * t18
      t72 = t17 * (-t20 + s16 * (s25 - t11)) + t14 * (-ddreal(9) * s16 +
     & t50)
      t73 = t2 * ((ddreal(16) * s26 + t28 + t31) * s16 - ddreal(6) * t61
     &) * t19
      t74 = s16 * t2 * t66 * (t10 * t2 + t7) * s34
      t42 = s16 * t14 * s25 - s16 * s56 * t42 + t72 * s26 + ((-(t16)**2 
     &* t30 - ddreal(9) * t14) * s16 + t3 * s16 * (s25 * t15 + (-t13 * t
     &56 - t16 * t62) * s26) + (t18 * ((-t52 + t59) * s34 + s16 * (s16 *
     & t65 - s25) + (ddreal(20) * s16 - ddreal(12) * s25 - ddreal(18) * 
     &s26) * s26) + t67 * (ddreal(9) * s25 + t63)) * s34 - ddreal(6) * t
     &13 * t57 * t16 - ddreal(6) * t2 * (-t4 * (-t24 + t16) + t29 - t36 
     &- t64 + t20 * t10) * t18) * s56 + t41 * ((-ddreal(16) * t15 - ddre
     &al(10) * t17) * s16 + t18 * (t51 * (-t39 + t52 - t7) + (-ddreal(45
     &) * s16 - ddreal(36) * s26) * s16 + ddreal(18) * t61 + ddreal(3) *
     & t68) + ddreal(3) * s34 * (t47 * (ddreal(16) * t2 + t7) + t19) + d
     &dreal(6) * s16 * ((-s16 * t57 - t69) * s26 + t44)) + t42 * ((s34 *
     & t10 + s25 - t11 - t54) * s16 + t71) + t43 * ((-t13 * t34 - ddreal
     &(10) * t16) * s16 + t3 * s16 * (-t69 + t24) - ddreal(6) * t19 + s3
     &4 * (s16 * (t7 + t63) + s34 * (-ddreal(12) * s16 + ddreal(12) * s2
     &6 + t70))) - t9 * (-s26 + t5) * (t18)**2 + t73 + t74
      t61 = t32 * t4 + t49
      t62 = t5 + t32
      t63 = t13 * s26
      t65 = s16 * t36
      t69 = ddreal(10) * s16
      t72 = t63 * t45
      t73 = ddreal(1) / s16
      t21 = ddreal(1) / (t21)**2
      t25 = ddreal(1) / (t25)**2
      t14 = (-t4 * t13 * (t16 * (s16 * t56 + s26 * t57) + t14) + t42 - t
     &71 * t66 * ((t59 + t7 + t53) * s16 - t58 - t60)) * MB1011(1) * t21
     & - t30 * t40 + ddreal(3) * t17 + ddreal(3) * t43 + ddreal(3) * t44
     & + (t37 + (t28 + t27) * s25 - t29) * s26 + (-(s26 * t26 + ddreal(1
     &0) * s25 + t35) * s34 + (s25 * t26 + s26 * t34 + t35) * s16 + t30 
     &* t33 - t16) * s34 + ((t39 + t31 + t38) * s56 + (s16 * t34 - t23 +
     & t31) * s34 + ddreal(9) * t16 - t26 * (t18 + t13) + t4 * (t37 + t3
     &6) + t28 * s25) * s56 + t5 * t36 + s25 * (t4 * t40 - (-(-s25 + t5)
     & * s26 - t13 * t3 + t68) * s26 - (-t4 * t61 * s34 + (t12 + t11 + t
     &50) * s16 + t4 * t58 + t64) * s34 - (-(t38 * t4 - s25) * s56 - (-t
     &28 - t7) * s34 - t3 * (-s26 * s34 + t13 + t18 + t20) + t60 + t68) 
     &* s56 - t65) * MB1110(1)
      t26 = ddreal(1) / t1
      t10 = (((-t15 * t45 - t17 * t45) * s16 + ((((s34 - ddreal(3) * t32
     & - t52) * s34 + (s25 * t10 + ddreal(12) * s26 + t69) * s16 + ddrea
     &l(3) * t58 + t70 * s26) * s34 - s25 * t36 - t13 * (t70 + t69) - t1
     &7 - ddreal(3) * ((s25 + t49) * s26 + (s25 + t49) * s25 + ddreal(6)
     & * t13) * s26 - ddreal(3) * t65) * s34 + t67 * (t52 + t23)) * s34 
     &+ (-t41 * (-t5 * s34 + t18 + t46) + ddreal(6) * (t18 * t62 + t47 *
     & (s26 + t5)) * s34 - ddreal(6) * t72 - ddreal(3) * (t13 * t45 + t1
     &6 * t45) * s16 - ddreal(3) * t18 * (t18 + (t49 + t23) * s25 + t48 
     &+ t58)) * s56 - ddreal(3) * t41 * ((s26 * t45 - s34 * (t49 + t23) 
     &+ t46) * s16 + t18 * (-s34 + t49 + t32)) - ddreal(3) * t72 * t2) *
     & MB1111D2(0) - s25 * s34 * (t4 * (((-t62 + s34) * s34 + t47) * s56
     & + t63) + (t16 + t13) * s16 + ((t61 - s34) * s34 - (t49 + t6) * s1
     &6 - t55 - t58) * s34 + t38 * t41) * MB1101(0)) * t25 * t26
      t10 = t73 * (t10 + t18 * MB1011(0) * t21)
      result = t25 * ((-t4 * (-s25 * t16 + t15 - t19) + (t16 - ddreal(3)
     & * t46 + t36) * s26 + (t49 * s25 - (t11 + t8 + t50) * s34 + t33 * 
     &t4 + t48) * s34 + ((t8 + t27 + s56) * s56 + (s25 + t51 + t6) * s25
     & + ddreal(3) * t16 - ddreal(3) * t18 - ddreal(3) * t46 + t11 * s34
     &) * s56 + t24 * t22) * I300s16(1) + t14 * t73 + t21 * ((t4 * (t18 
     &* (s16 + t23) + t15) + (s16 * (t52 - t27) + (s26 - s25 + t53) * s2
     &6) * s26 + ((-s34 * t30 - s25 + s56 + t53 + t8) * s56 + s16 * (-dd
     &real(9) * s34 + t52) + (-ddreal(10) * s34 + t28 + t8) * s26 + t3 *
     & s34 * (-s25 + s34) - t55) * s56 - t44 - t2 * (t22 * t3 + t54) * s
     &34) * MB1001(0) + (t9 * s34 - t2 * (t7 + t6 + t5)) * s34 + ((-s34 
     &* t3 + s56 + ddreal(3) * t2) * s56 - (t12 + t11 + t7) * s34 + ddre
     &al(3) * t18 + ddreal(3) * t16 + ddreal(3) * t13 + t11 * s26) * s56
     & + t15 + t17 + ddreal(3) * t20 * t2) * t1) / ddreal(4) + t10 / ddr
     &eal(2)

           intHs16s25s26s34s56x1132D6eps0 = result
       end function intHs16s25s26s34s56x1132D6eps0

       function intHs16s25s26s34s56x1132D6eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1132D6eps1
           type(dd_complex) ::  t1,t2

           type(dd_complex) :: result

      t1 = s16 + s26 + s56
      t2 = ddreal(1) / (t1)**2
      result = ddreal(1) / s16 * t2 * (s34 / ddreal(2) + t1 / ddreal(4) 
     &+ (s34)**2 * MB1011(1) / ddreal(2))

           intHs16s25s26s34s56x1132D6eps1 = result
       end function intHs16s25s26s34s56x1132D6eps1

       function intHs16s25s26s34s56x1141D6eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1141D6eps0
           type(dd_complex) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           type(dd_complex) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           type(dd_complex) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           type(dd_complex) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           type(dd_complex) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           type(dd_complex) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           type(dd_complex) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           type(dd_complex) ::  t175,t18,t19,t2,t20,t21,t22,t23,t24,t25,t26,t27
           type(dd_complex) ::  t28,t29,t3,t30,t31,t32,t33,t34,t35,t36,t37,t38
           type(dd_complex) ::  t39,t4,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49
           type(dd_complex) ::  t5,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t6
           type(dd_complex) ::  t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t7,t70
           type(dd_complex) ::  t71,t72,t73,t74,t75,t76,t77,t78,t79,t8,t80,t81
           type(dd_complex) ::  t82,t83,t84,t85,t86,t87,t88,t89,t9,t90,t91,t92
           type(dd_complex) ::  t93,t94,t95,t96,t97,t98,t99

           type(dd_complex) :: result

      t1 = ddreal(3)
      t2 = (s16)**2
      t3 = (t2)**2
      t4 = s16 * t2
      t5 = t2 * t3
      t6 = s16 * t3
      t7 = t1 * s25
      t8 = ddreal(19) * s16
      t9 = ddreal(25)
      t10 = ddreal(17) * s16
      t11 = t9 * s25
      t12 = s16 - s25
      t13 = s26 * t12
      t14 = ddreal(5)
      t15 = ddreal(4)
      t16 = (s26)**2
      t17 = (t16)**2
      t18 = s26 * t17
      t19 = s26 * t16
      t20 = t1 * s16
      t21 = t14 * s25
      t22 = ddreal(15) * t2
      t23 = ddreal(2)
      t24 = t23 * s26
      t25 = ddreal(10) * s25
      t26 = ddreal(15) * s16
      t27 = s16 - s34
      t28 = s26 * t27
      t29 = (s34)**2
      t30 = (t29)**2
      t31 = s34 * t30
      t32 = s34 * t29
      t33 = t2 + t29
      t34 = t23 * s34 * t12
      t35 = ddreal(42)
      t36 = -t32 + t4
      t37 = (s56)**2
      t38 = (t37)**2
      t39 = s56 * t37
      t40 = s56 * t38
      t41 = t19 + t39
      t42 = ddreal(6) * s25
      t43 = ddreal(13) * s16
      t44 = s16 * s25
      t45 = s16 + s25
      t46 = (-s34 + t45 + s26) * s26 + (t24 + s16 + s25 - s34 + s56) * s
     &56 + t44
      t47 = ddreal(12) * s25
      t48 = t14 * t2
      t49 = ddreal(9)
      t50 = ddreal(6) * s16
      t51 = t49 * s25
      t52 = s16 * t49
      t53 = (-t52 - t42) * s25
      t54 = ddreal(12) * s16
      t55 = t7 + t54
      t56 = s25 * t55 - t48
      t57 = t7 - s16
      t58 = s16 + s26
      t59 = t23 * s25
      t60 = ddreal(13) * t2
      t61 = ddreal(11) * s16
      t62 = (s25)**2
      t63 = s25 * t62
      t64 = t1 * t19
      t65 = ddreal(35)
      t66 = ddreal(14) * s25
      t67 = ddreal(34) * s16
      t68 = ddreal(20) * s25
      t69 = ddreal(40) * s16
      t70 = (t69 + t51) * s25
      t71 = ddreal(38) * t2
      t72 = ddreal(46) * s16
      t73 = ddreal(7)
      t74 = t73 * s25
      t75 = ddreal(10) * t2
      t76 = ddreal(12) * s26
      t77 = t1 * s26
      t78 = t77 + s16
      t79 = t15 * t56
      t80 = -t57
      t81 = ddreal(20) * s16
      t82 = t15 * s25
      t83 = ddreal(53) * s16
      t84 = t15 * s16
      t85 = t84 + t59
      t86 = t85 * s25
      t87 = ddreal(15) * t19
      t88 = (ddreal(55) * s16 + t66) * s25
      t89 = ddreal(66) * t2
      t90 = ddreal(138) * s16
      t91 = ddreal(60) * s25
      t92 = t9 * s16
      t93 = ddreal(50) * t2
      t94 = t35 * t4
      t95 = ddreal(54) * s26
      t96 = ddreal(24)
      t97 = t96 * s26
      t98 = t3 * s25
      t99 = s26 * (t48 + t53)
      t100 = ddreal(6) * t99
      t101 = s16 * ((-t51 - t50) * s25 + t48)
      t102 = ddreal(12) * t62
      t103 = t84 + s25
      t104 = t103 * s25
      t105 = ddreal(19) * t2
      t106 = t14 * s16
      t107 = s25 + t106
      t108 = ddreal(18) * t2
      t109 = t4 * s25
      t110 = s16 * t62
      t111 = ddreal(15) * s26
      t112 = ddreal(23) * s16
      t113 = t2 * s25
      t114 = t80 * t16
      t115 = t1 * t62
      t116 = s25 + s34
      t117 = (t58)**2
      t105 = ((-(t1 * t116 - s16) * s56 - (s25 + t111 + t10) * s34 + t14
     & * (s26 * t80 + t2) + ddreal(12) * t29 - ddreal(12) * t44 - t115) 
     &* s56 - s26 * t79 - t23 * s34 * (-(t97 + t112 + t25) * s34 + (t59 
     &+ t67 + t111) * s26 + t104 + t105) - ddreal(18) * t113 + ddreal(10
     &) * t114 - ddreal(18) * t32 + ddreal(10) * t4 - t102 * s16) * s56 
     &- (ddreal(18) * t110 - t100) * s16 + ddreal(10) * t19 * t80 - ddre
     &al(12) * t109 + ddreal(10) * t3 + ddreal(12) * t30 - (((t72 + t95 
     &+ t11) * s34 - (ddreal(72) * s26 + t90 + t91) * s26 - t88 - t89) *
     & s34 + (s25 * t107 + t108) * s25 + ddreal(30) * t19 + t94 + ddreal
     &(6) * s26 * ((s25 + t10) * s26 + t104 + t105)) * s34 - ddreal(6) *
     & t16 * t56
      t79 = t105 * s56 - s16 * ((-t101 * t15 - t100) * s26 + t102 * t2) 
     &- t1 * (t31 + t98) + t14 * (t17 * t80 + t6) - (((-(t97 + t51 + t10
     &) * s34 + (ddreal(92) * s16 + ddreal(50) * s25 + t95) * s26 + t70 
     &+ t71) * s34 - ((t7 + t92) * s25 + t93) * s25 - ((ddreal(48) * s26
     & + t90 + t91) * s26 + t23 * (t89 + t88)) * s26 - t94) * s34 + t58 
     &* ((ddreal(16) * t2 + t86) * s25 + ddreal(23) * t4 + s26 * ((t82 +
     & t83) * s26 + (t81 + t42) * s25 + ddreal(61) * t2) + t87)) * s34 -
     & t19 * t79
      t88 = s16 * t16
      t89 = s16 + s26 + s56
      t90 = s16 + s25 + s26 - s34 + s56
      t91 = s16 + s26 - s34 + s56
      t94 = ddreal(16) * s16
      t95 = ddreal(11) * s25
      t97 = -ddreal(22)
      t100 = ddreal(18) * s26
      t105 = ddreal(10) * s16
      t118 = -s25 * t23 * t45 + t2
      t86 = t2 - t86
      t119 = t59 + s16
      t120 = t23 * s16
      t121 = t59 + t106
      t122 = -t74 + t105
      t123 = ddreal(13) * s25
      t124 = ddreal(10) * s26
      t125 = ddreal(20) * t4
      t126 = ddreal(60) * t2
      t127 = ddreal(29) * s25
      t128 = ddreal(60) * s16
      t129 = t77 + t120
      t130 = ddreal(40) * t2
      t131 = ddreal(20) * t2
      t83 = (t123 + t83) * s25
      t132 = ddreal(90) * t2
      t133 = ddreal(8) * s25
      t134 = ddreal(70) * t4
      t135 = ddreal(26) * s16
      t136 = ddreal(30) * s26
      t137 = s16 * s26
      t138 = ddreal(6) * s26
      t139 = s25 * t107
      t140 = s16 * t58
      t141 = t73 * s26
      t142 = ddreal(40) * s25
      t143 = s25 + s26
      t104 = ((t37 + (-s34 * t73 + t136 + t21 + t92) * s25 + ddreal(15) 
     &* s26 * t58 + ddreal(10) * s34 * t27) * s56 + (t15 * ((-t84 - s25 
     &- t141) * s25 + ddreal(10) * t140) - (ddreal(20) * s26 + t133 + t9
     &2) * t23 * s34) * s34 + t2 * (t142 - t105) + ddreal(20) * (t139 + 
     &t16) * s26 + ddreal(20) * t110 + ddreal(20) * t32 + ddreal(30) * t
     &119 * t16) * s56 + ddreal(30) * (-s16 * t86 + (s26 * t119 + t139) 
     &* s26) * s26 + ddreal(30) * t113 * t45 - ddreal(15) * t3 + ddreal(
     &15) * t17 - ddreal(15) * t30 + (((ddreal(60) * t58 + t127) * s34 -
     & (ddreal(150) * s16 + ddreal(48) * s25 + ddreal(60) * s26) * s26 -
     & t132 - t83) * s34 + s16 * ((-t95 - t50) * s25 + t126) + ddreal(6)
     & * t122 * t16 - t63 + t76 * (t75 - t104)) * s34
      t52 = t104 * s56 + t15 * t31 + t17 * (ddreal(15) * t119 + t138) + 
     &t23 * s34 * (t29 * ((ddreal(29) * s26 + t42 + t135) * s25 + ddreal
     &(30) * t16 + ddreal(30) * t2 + t128 * s26) + t58 * (s16 * ((t84 - 
     &t21) * s25 + t131) - (-t23 * t122 * s26 - (-t42 - t105) * s25 - t1
     &30) * s26 - t63)) + t3 * (-t52 + t25) + ddreal(20) * s25 * (t107 *
     & t19 + t109) - t29 * (t29 * (t136 + t92 + t47) + ((t82 + t92) * s2
     &5 + ddreal(58) * t2) * s25 + ((ddreal(40) * s26 + ddreal(6) * t133
     & + ddreal(6) * t92) * s26 + t23 * (t132 + t83)) * s26 + t134) - dd
     &real(30) * t137 * (s16 * t118 + s26 * t86)
      t52 = t1 * (t119 * t18 + t40 * (t143 * t23 + s16)) + t14 * s25 * (
     &t107 * t17 + t98) - t23 * (-t31 * (t24 + s16) + t5) + s25 * t6 + (
     &t18 + t4 * ((t105 + t68) * s25 - t2 * t49)) * s26 + s34 * (((-t129
     & * (t14 * t58 + t82) * s34 + ddreal(20) * t19 + ddreal(20) * t4 + 
     &s26 * ((t128 + t127) * s26 + (ddreal(52) * s16 + t47) * s25 + t126
     &) + t44 * (t112 + t25)) * s34 - t58 * (((t82 + t54) * s25 + ddreal
     &(21) * t2) * s25 + ((ddreal(8) * t121 + t124) * s26 + (ddreal(37) 
     &* s16 + t123) * s25 + t93) * s26 + t125)) * s34 + t117 * (s16 * ((
     &-t7 + t106) * s25 + t75) + t122 * t16 - t63 + t24 * (t120 - s25) *
     & t121)) + t52 * s56 - t88 * (t118 * t26 + t124 * t86)
      t83 = s16 * t12
      t86 = t120 - t21
      t98 = t50 + t21
      t104 = ddreal(14) * s16
      t107 = (-t104 - t25) * s25 + t48
      t109 = (-t43 - t21) * s25 + t2 * t23
      t110 = ddreal(31) * s16
      t112 = (t110 + t21) * s25
      t118 = t1 * t2 + t112
      t119 = t7 + t120
      t121 = ddreal(8) * s16
      t122 = ddreal(6) * t2
      t126 = t96 * t4
      t127 = t84 - s25
      t128 = (t20 - t59) * t127
      t132 = t128 * s26
      t139 = s16 * t63
      t144 = t14 * t4
      t145 = ddreal(75)
      t146 = ddreal(48) * t2
      t147 = t20 + s25
      t148 = -s25 + t26
      t149 = ddreal(30) * s16
      t150 = ddreal(21) * s25
      t151 = (-t42 + t106) * s25
      t152 = ddreal(37) * t2
      t153 = t152 + t151
      t154 = ddreal(27)
      t155 = t154 * s16
      t156 = t155 - t11
      t157 = t15 * t4
      t158 = ddreal(60) * t4
      t159 = t132 * t14
      t160 = ddreal(65) * s16
      t161 = ((-t82 - t20) * s25 + t146) * s25
      t162 = ddreal(79) * t4
      t163 = (-ddreal(78) * s16 + ddreal(63) * s25) * s25 - ddreal(453) 
     &* t2
      t164 = ddreal(56) * s25
      t165 = -ddreal(124) * s16
      t166 = t164 + t165
      t167 = ddreal(36) * s26
      t168 = ddreal(12) * t119
      t169 = t118 * t14
      t170 = s16 * t109
      t171 = t22 * t107
      t172 = t9 * t3
      t55 = t1 * (s16 * t55 + (-s26 * t49 + t105 - t133) * s26) * t31
      t55 = s16 * ((((((-t168 - t141) * s26 - t169) * s26 + ddreal(20) *
     & t170) * s26 + t171) * s26 + t157 * t86 * t98) * s26 - t172 * t62)
     & + ddreal(11) * t5 * t12 - s34 * ((((-t29 * (-t138 + t106) - s16 *
     & ((s16 * t145 - t7) * s25 + ddreal(105) * t2) - ((t1 * t156 - t167
     &) * s26 + ddreal(6) * t153) * s26) * s34 - s16 * (((-ddreal(36) * 
     &s16 + t123) * s25 - ddreal(154) * t2) * s25 - ddreal(160) * t4) - 
     &(-ddreal(6) * t162 - ddreal(6) * t161 + t87) * s26 - t16 * (s26 * 
     &t166 + t163)) * s34 - t58 * ((((t20 - t42) * s25 + t71) * s25 + dd
     &real(108) * t4) * s25 + ddreal(135) * t3 + ((t14 * t148 * s26 + (-
     &t82 + t160) * s25 + ddreal(285) * t2) * s26 + ((t67 - t51) * s25 +
     & ddreal(178) * t2) * s25 + ddreal(345) * t4) * s26)) * s34 + t117 
     &* ((((s16 + t42) * s25 - ddreal(29) * t2) * s25 + t4 * t49) * s25 
     &+ (t23 * (((-t104 + t133) * s25 - ddreal(23) * t2) * s25 + t158) +
     & t159) * s26 + ddreal(60) * t3)) - t55
      t87 = -ddreal(43)
      t104 = t128 * t124
      t123 = ddreal(120) * t3
      t124 = ddreal(90) * s16
      t173 = (t124 - t42) * s25 + ddreal(540) * t2
      t174 = ddreal(10) * t148 * s26
      t175 = ddreal(10) * t19
      t18 = s16 * (s26 * t171 - t118 * t175 - t142 * t3 - ddreal(50) * t
     &4 * t62 + t6 * t96 - ddreal(21) * t18) - t1 * ((t29 - (t155 - t11 
     &- t100) * s26 - t151 - t152) * s34 + t162 - t175 + t161 + s26 * (-
     &s26 * t166 / ddreal(2) - t163 / ddreal(3))) * t32 + ddreal(30) * t
     &88 * (-t119 * t16 + t170) - ((-t32 * (s26 * t154 - t26 + t47) + t1
     & * t63 * t45 - t2 * ((ddreal(143) * s16 + ddreal(36) * s25) * s25 
     &+ ddreal(240) * t2) - ((t174 + t173) * s26 + ((t149 - t51) * s25 +
     & ddreal(243) * t2) * s25 + ddreal(630) * t4) * s26) * s34 + t58 * 
     &(((t2 * t87 - t53) * s25 - ddreal(14) * t4) * s25 + (t15 * (((-t12
     &1 + t42) * s25 - ddreal(31) * t2) * s25 + t158) + t104) * s26 + t1
     &23)) * s34
      t53 = (t82 - t120) * s25
      t87 = ddreal(39) * t2
      t110 = s16 * (s16 * ((-t160 - t11) * s25 + t75) + ((-s26 * t65 - d
     &dreal(30) * t119) * s26 - t169) * s26) - t30 * t49 + (((-t110 + t1
     &11 + t66) * s34 + s16 * (s26 * t145 + t124) + (-s26 * t14 - s25 + 
     &t26) * s25) * s34 + (t87 - t53) * s25 - t158 - t159) * s34
      t111 = s16 * ((ddreal(21) * s26 + t168) * s26 + t112) + t1 * t36 +
     & s34 * (-s34 * t148 + t128)
      t112 = s25 * t5 + (t17 * (t16 + t118) - t3 * (-t62 * t9 + ddreal(1
     &1) * t83)) * s26
      t69 = t58 * ((t15 * (t1 * t62 * t12 + t2 * (ddreal(26) * s25 + t14
     &9)) - t64) * s26 + t2 * ((s25 * t96 + ddreal(50) * s16) * s25 + t1
     &30) - t16 * (-(t67 - t66) * s26 - (-t150 + t69) * s25 - ddreal(117
     &) * t2) - t139) * t32
      t54 = t117 * (t1 * t63 * t12 + ddreal(30) * t4 * t45 + ((t148 * s2
     &6 + (t81 - s25) * t147) * s26 + ((-t7 + t54) * s25 + t146) * s25 +
     & t145 * t4) * s26 + t60 * t62) * t29
      t81 = t58 * t117 * (((-t2 * t73 + t23 * t62) * s25 + t144) * s25 +
     & ddreal(12) * t3 - (-t132 - ((t82 - t121) * s25 - t122) * s25 - t1
     &26) * s26 - t139) * s34
      t18 = s16 * t39 * t38 + s16 * t112 - t55 * s56 - t30 * (t1 * s26 *
     & (s16 * ((-s25 + t92) * s25 + t2 * t65) + t153 * s26) + t156 * t19
     & + t4 * (s25 * t65 + t149) + t49 * (t2 * t62 - t17)) - t37 * t18 -
     & t38 * t110 - t39 * (s16 * (-s25 * t134 + ddreal(20) * t109 * t137
     & - ddreal(10) * t118 * t16 - ddreal(40) * t119 * t19 - t17 * t65 -
     & t62 * t93 + t172) + t31 * t49 + ((((t155 - t167 - t11) * s34 + (t
     &150 - t135) * s25 + (t164 + t165 + t136) * s26 - ddreal(151) * t2)
     & * s34 + ((-t7 + t105) * s25 + ddreal(81) * t2) * s25 + ddreal(210
     &) * t4 + s26 * (t174 + ddreal(2) / ddreal(3) * t173)) * s34 + (((-
     &t59 - t61) * s25 + t2 * t9) * s25 + ddreal(46) * t4) * s25 + (-t15
     & * ((-t87 + t53) * s25 + t158) - t104) * s26 - t123) * s34) + t40 
     &* t111 - t78 * (t120 - s26) * t29 * t30 + t69 + s16 * (t84 + t42 +
     & t141) * t37 * t38 - t54 + t81
      t38 = t16 + t62
      t39 = ddreal(6) * t140
      t40 = t120 + t143
      t49 = t2 + t62
      t53 = t20 + t143 - s34
      t54 = (t12)**2
      t55 = t54 * s26
      t69 = s16 * t49
      t81 = t23 * s16 * t116
      t87 = t40 * t32
      t49 = t2 * t49
      t54 = t54 * t16
      t93 = t63 - t4
      t104 = s16 * t93
      t93 = t1 * t93 + t44 * (t82 + t106)
      t110 = t120 + s25
      t111 = t14 * t44 * t12
      t112 = t1 * t63
      t108 = s25 * (-t84 - t42) + t108
      t105 = t133 + t105
      t114 = t114 * t110
      t116 = t114 - s25 * (s25 * (t75 + t115) + t144) - ddreal(6) * t104
      t80 = s25 * (s26 * t116 + s34 * (s34 * (-s34 * (t96 * s16 * t143 +
     & s26 * (ddreal(15) * s25 + t138) + t102 + t131) + s25 * (s25 * (t1
     &33 + t26) + t2 * t96) + s26 * (s26 * (t24 + ddreal(6) * t147) + s2
     &5 * (t51 + t92) + ddreal(36) * t2) + t125) + t23 * t62 * (s25 * (s
     &16 - s25) + t2) - t4 * t105 - (s26 * (s26 * t127 + t108) + t111 - 
     &t112 + t126) * s26) + s56 * (s56 * (-s56 * (s25 * (t7 - s34 + t106
     &) - t23 * t33 + t84 * s34) - s34 * (s16 * (ddreal(18) * s16 + t76)
     & + s25 * (-t84 - t77 - t42)) + ddreal(6) * t29 * t53 + ddreal(6) *
     & t4 - ddreal(6) * t63 + t77 * t80 * t110 - t44 * t105) - t1 * (t32
     & * (s26 * t15 + t121 + t21) + (t62)**2 - t114) - ddreal(6) * t104 
     &+ ddreal(6) * t30 - t15 * t93 * s26 - s34 * (-s34 * (s25 * (t76 + 
     &t51 + t92) + ddreal(36) * t140 + ddreal(6) * t16) + s26 * (s26 * t
     &1 * t127 + t15 * t108 / ddreal(2)) + t111 - t112 + t126) - t113 * 
     &(t106 + t25)))
      t92 = ddreal(1) / s16
      t90 = ddreal(1) / (t90)**2
      t46 = ddreal(1) / (t46)**2
      t89 = ddreal(1) / (t89)**2
      t1 = (-t14 * t2 * (t19 * (-s16 * t107 - s26 * t109) + t3 * t62) + 
     &t23 * s16 * (t16 * (t4 * t86 * t98 - t119 * t17) + t4 * t3) - t18 
     &- t1 * (t1 * ((s16 * t103 - t16) * s26 + t113) + t16 * (-t82 + t10
     &6) + t157) * t31) * MB1011(1) * t92
      t1 = t89 * (t52 * MB1001(0) - t23 * t88 * (t101 + t99) - (t19 * (-
     &s26 * t57 - t56) - t4 * ((t47 + t20) * s25 - t48)) * s26 - s34 * (
     &(((-t78 * s34 + t14 * s16 * t45 + (t76 + t51 + t10) * s26) * s34 -
     & s16 * ((t74 + t26) * s25 + t75) - ddreal(18) * t19 - s26 * ((t72 
     &+ t11) * s26 + t70 + t71)) * s34 + t58 * (((t7 + t61) * s25 + t22)
     & * s25 + ddreal(12) * t19 + ddreal(10) * t4 + s26 * ((t67 + t68) *
     & s26 + (s16 * t65 + t66) * s25 + ddreal(32) * t2))) * s34 - t117 *
     & (t14 * t2 * t45 + ((s25 + t61) * s26 + (t59 + t50) * s25 + t60) *
     & s26 + t45 * t62 + t64)) - t79 * s56 - t5 + t115 * t3 + t1) + (-t2
     &3 * s25 * (-t16 * t93 + t30 * (t82 + t77 + t106 - s34) + t6) - t80
     & + t139 * (t48 + t115)) * MB1110(1) * t92
      t3 = ddreal(1) / t91
      t3 = t3 * (t92 * (-t32 * MB1011(0) * t89 + (-s25 * s34 * (-t23 * (
     &s16 * (s34 * t58 * (t120 - s25 + s26) + t113 - t55) - s56 * (-t113
     & * t23 + (-s16 * (t20 + t24 - s25) + t53 * s34) * s34 + t55 + t69)
     & + t87) + t29 * (s25 * t23 * t58 + t29 + t38 + t39) + t37 * (-t81 
     &+ t2 + t62 + t29) + t49 + t54) * MB1101(0) + (-s26 * t12 + s34 * (
     &t40 - s34) - (s16 - s25 - s34) * s56 - t83) * (-t23 * (s16 * (t113
     & - t55) + t87) + s34 * (s34 * (s25 * t129 + t29 + t38 + t39) - t58
     & * (s16 * (t84 - t59) + s26 * t110 + t62)) + s56 * (s56 * (-t81 - 
     &s34 * (s25 - s34) + t2 + t62) + t23 * (-t32 + t55 + t69) - s34 * (
     &-s25 * t12 + s26 * t85 - s34 * (t7 + t24 + t50) + t122) - t113 * t
     &15) + t49 + t54) * MB1111D2(0)) * t3 * t46 - t63 * MB1110(0) * t90
     &) + I300s16(0))
      result = t46 * ((t14 * t36 + t23 * (((-t7 - t8) * s25 + t2) * s26 
     &+ ((-s26 * t9 - t7 - t8) * s25 - ddreal(21) * t16 - ddreal(17) * t
     &28 + t33 - t34) * s56) - ddreal(14) * t41 + ((t24 + t26 + t25) * s
     &34 + (t20 - t21) * s25 - t13 * t15 + ddreal(17) * t16 - t22) * s34
     & + (-t11 - t10) * t16 + (-s26 * t35 - t11 - ddreal(17) * t27) * t3
     &7 - t44 * (t43 + t42)) * t92 + (t23 * t36 - ddreal(6) * t41 + ((-s
     &16 * t73 - t95) * s26 + (-t7 - t94) * s25 + t2) * s26 + ((t82 + s2
     &6 + t50) * s34 + (s25 - t50) * s16 + t16 * t73 - t23 * (t62 + t13)
     &) * s34 + ((-t27 * t73 - t100 - t95) * s56 + (s26 * t97 - t7 - t94
     &) * s25 - ddreal(18) * t16 - ddreal(14) * t28 + t33 - t34) * s56 -
     & t44 * (t7 + t106)) * I300s16(1) + t1 * t90) / ddreal(12) - t3 / d
     &dreal(6)

           intHs16s25s26s34s56x1141D6eps0 = result
       end function intHs16s25s26s34s56x1141D6eps0

       function intHs16s25s26s34s56x1141D6eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1141D6eps1
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t2,t3,t4,t5,t6,t7
           type(dd_complex) ::  t8,t9

           type(dd_complex) :: result

      t1 = (s16)**2
      t2 = ddreal(3) * s25
      t3 = ddreal(11) * s25
      t4 = (s26)**2
      t5 = t4 + t1
      t6 = s16 + s26
      t7 = (s25)**2
      t8 = ddreal(6)
      t9 = (s34)**3
      t10 = s16 + s26 + s56
      t11 = s16 + s25 + s26 - s34 + s56
      t12 = s16 + s26 - s34 + s56
      t13 = ddreal(1) / s16
      t12 = ddreal(1) / t12
      t11 = ddreal(1) / (t11)**2
      t10 = ddreal(1) / (t10)**2
      t14 = -ddreal(1) / ddreal(12)
      result = t12 * ((s25 * t7 * MB1110(1) * t11 + t9 * MB1011(1) * t10
     &) * t13 - I300s16(1)) / ddreal(6) + t14 * (t8 * ((s56)**3 + s16 * 
     &t1 + s26 * t4) + ((ddreal(18) * s16 + t3) * s26 + (ddreal(22) * s1
     &6 + t2) * s25 + ddreal(18) * t1) * s26 + (-(ddreal(4) * s25 + t6) 
     &* s34 - ddreal(7) * t5 + ddreal(2) * (-s16 - s26 + s25) * s25 - dd
     &real(14) * s16 * s26) * s34 + ((-ddreal(7) * s34 + t3 + ddreal(18)
     & * t6) * s56 + ddreal(22) * s25 * t6 - ddreal(14) * s34 * t6 + (-d
     &dreal(2) * s25 - s34) * s34 + ddreal(18) * t5 + ddreal(3) * t7 + d
     &dreal(36) * s16 * s26) * s56 + ddreal(2) * t9 + s16 * s25 * (ddrea
     &l(11) * s16 + t2)) * t13 * t10 * t11

           intHs16s25s26s34s56x1141D6eps1 = result
       end function intHs16s25s26s34s56x1141D6eps1

       function intHs16s25s26s34s56x1210D2eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1210D2eps0
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s16 + s25 + s26 - s34 + s56
      t1 = ddreal(1) / t1
      result = ((s16 + s26 + s56) * MB1011(1) + (s16 + s26 - s34 + s56) 
     &* MB1110(0) + (s16 - s25 + s26 - s34 + s56) * MB1110(1) - MB1001(0
     &)) * t1

           intHs16s25s26s34s56x1210D2eps0 = result
       end function intHs16s25s26s34s56x1210D2eps0

       function intHs16s25s26s34s56x1210D2eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1210D2eps1
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s16 + s25 + s26 - s34 + s56
      t1 = ddreal(1) / t1
      result = t1 * ((s16 + s26 - s34 + s56) * MB1110(1) - ddreal(1))

           intHs16s25s26s34s56x1210D2eps1 = result
       end function intHs16s25s26s34s56x1210D2eps1

       function intHs16s25s26s34s56x1211D2eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1211D2eps0
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s34 - s16 - s26 - s56
      t1 = ddreal(1) / t1
      result = -ddreal(1) / s16 * ((s34 * MB1101(0) - (s16 + s26 + s56) 
     &* MB1111D2(0)) * t1 - MB1110(0))

           intHs16s25s26s34s56x1211D2eps0 = result
       end function intHs16s25s26s34s56x1211D2eps0

       function intHs16s25s26s34s56x1211D2eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1211D2eps1

           type(dd_complex) :: result

      result = MB1110(1) / s16

           intHs16s25s26s34s56x1211D2eps1 = result
       end function intHs16s25s26s34s56x1211D2eps1

       function intHs16s25s26s34s56x1211D4eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1211D4eps0
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           type(dd_complex) ::  t20,t3,t4,t5,t6,t7,t8,t9

           type(dd_complex) :: result

      t1 = s16 - s25
      t2 = s16 + s25
      t3 = (s16)**2
      t4 = t2 * s26
      t5 = ddreal(2) * t3
      t6 = (s25)**2
      t7 = s16 * s25
      t8 = (s34)**2
      t9 = t8 + t6 + t3
      t10 = (t1)**2 * s26
      t11 = t5 * s25
      t5 = -t11 + t9 * s16 - (t5 + t4) * s34 + (-s16 * s34 - s25 * s34 +
     & t3 + t6 - ddreal(2) * t7) * s56 + t10
      t12 = -ddreal(2) * s34 * t2 - ddreal(2) * t7 + t9
      t13 = ddreal(2) * s26
      t14 = (-s34 + t2 + s26) * s26 + (t13 + s16 + s25 - s34 + s56) * s5
     &6 + t7
      t15 = s16 + s26 + s56
      t16 = ddreal(3) * s16
      t17 = (t6 + t3) * s16
      t18 = t1 * s16
      t19 = ddreal(1) / t14
      t12 = ddreal(1) / t12
      t14 = ddreal(1) / t14
      t20 = -ddreal(1) / ddreal(2)
      result = t20 * (-t15 * (s16 + s26 - s34 + s56) * MB1111D2(0) * t14
     & + (((ddreal(2) * s34 * (s56)**2 - ddreal(2) * s25 * t3 + (ddreal(
     &2) * (s26)**2 + (s16 - t13) * s34 - ddreal(2) * t18 + t4) * s34 + 
     &(-ddreal(2) * t7 - ddreal(2) * t8 + (ddreal(4) * s26 + t2) * s34 +
     & t3 + t6) * s56 + t10 + t17) * MB1110(1) + s34 * (-s26 * t1 + (s16
     & - s26) * s34 - (s16 - s25 + s34) * s56 - t18) * MB1101(0)) * s25 
     &+ t5 * (-s16 * I300s16(1) + MB1001(0)) - t15 * (-t11 + (-s16 * (-s
     &25 + t16 + t13) + (-s34 + t16 + s25 + s26) * s34) * s34 + (t9 - dd
     &real(2) * s16 * (s25 + s34)) * s56 + t10 + t17) * MB1011(1)) * t12
     & * t19) + t5 * t12 * t19

           intHs16s25s26s34s56x1211D4eps0 = result
       end function intHs16s25s26s34s56x1211D4eps0

       function intHs16s25s26s34s56x1211D4eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1211D4eps1

           type(dd_complex) :: result

      result = ddreal(0)

           intHs16s25s26s34s56x1211D4eps1 = result
       end function intHs16s25s26s34s56x1211D4eps1

       function intHs16s25s26s34s56x1212D4eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1212D4eps0
           type(dd_complex) ::  t1,t10,t11,t12,t2,t3,t4,t5,t6,t7,t8,t9

           type(dd_complex) :: result

      t1 = s16 - s25
      t2 = ddreal(2) * s16
      t3 = t1 * s16
      t1 = t1 * s26
      t4 = (t2 + s25 + s26 - s34) * s34 - (s16 - s25 - s34) * s56 - t3 -
     & t1
      t5 = s16 + s25
      t6 = s16 * s25
      t7 = (s34)**2
      t8 = (s16)**2
      t9 = (s25)**2 - ddreal(2) * s34 * t5 - ddreal(2) * t6 + t7 + t8
      t10 = s16 + s25 - s34
      t11 = (-s34 + t5 + s26) * s26 + (ddreal(2) * s26 + t10 + s56) * s5
     &6 + t6
      t12 = s25 * s34
      t11 = ddreal(1) / t11
      t9 = ddreal(1) / t9
      result = ddreal(2) * t4 * t9 * t11 - t11 * ((s16 + s26 - s34 + s56
     &) * MB1111D2(0) + t9 * (t4 * (-s16 * I300s16(1) + MB1001(0)) + (s2
     &5 * (s16 + s25 + s26 - s34 + s56) * ((s16 - s26) * s34 - (s16 - s2
     &5 + s34) * s56 - t1 - t3) * MB1110(1) - (s16 + s26 + s56) * ((-t3 
     &- t1) * s16 + ((-ddreal(3) * s16 - ddreal(2) * s25 - s26 + s34) * 
     &s34 + ddreal(3) * t8 + s25 * t5 + (s25 + t2) * s26) * s34 - (-t2 *
     & s34 - t12 - t6 + t7 + t8) * s56) * MB1011(1)) / s16 + t12 * (t10 
     &+ ddreal(2) * s26 + ddreal(2) * s56) * MB1101(0)))

           intHs16s25s26s34s56x1212D4eps0 = result
       end function intHs16s25s26s34s56x1212D4eps0

       function intHs16s25s26s34s56x1212D4eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1212D4eps1

           type(dd_complex) :: result

      result = ddreal(0)

           intHs16s25s26s34s56x1212D4eps1 = result
       end function intHs16s25s26s34s56x1212D4eps1

       function intHs16s25s26s34s56x1213D6eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1213D6eps0
           type(dd_complex) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           type(dd_complex) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           type(dd_complex) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           type(dd_complex) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           type(dd_complex) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           type(dd_complex) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           type(dd_complex) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           type(dd_complex) ::  t175,t176,t177,t178,t18,t19,t2,t20,t21,t22,t23,t24
           type(dd_complex) ::  t25,t26,t27,t28,t29,t3,t30,t31,t32,t33,t34,t35
           type(dd_complex) ::  t36,t37,t38,t39,t4,t40,t41,t42,t43,t44,t45,t46
           type(dd_complex) ::  t47,t48,t49,t5,t50,t51,t52,t53,t54,t55,t56,t57
           type(dd_complex) ::  t58,t59,t6,t60,t61,t62,t63,t64,t65,t66,t67,t68
           type(dd_complex) ::  t69,t7,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79
           type(dd_complex) ::  t8,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t9
           type(dd_complex) ::  t90,t91,t92,t93,t94,t95,t96,t97,t98,t99

           type(dd_complex) :: result

      t1 = s16 - s25
      t2 = ddreal(2)
      t3 = (s25)**2
      t4 = (t3)**2
      t5 = s25 * t4
      t6 = s25 * t3
      t7 = t2 * t3
      t8 = ddreal(3) * s16 * t1 - t7
      t9 = ddreal(3) * s16
      t10 = -s25 - t9
      t11 = (s16)**2
      t12 = (t11)**2
      t13 = s16 * t11
      t14 = s16 * t12
      t15 = t2 * t11
      t16 = ddreal(3) * s25
      t17 = s16 + s26
      t18 = -t3 + t11
      t19 = s16 * s25
      t20 = ddreal(8)
      t21 = ddreal(7) * t11
      t22 = ddreal(10) * t12
      t23 = t20 * s16
      t24 = t2 * t13
      t25 = ((-s25 + t23) * s25 + t21) * s25 + t24
      t26 = t11 + t3
      t27 = t26 * t3
      t28 = (s26)**2
      t29 = s26 * t28
      t30 = ddreal(24) * s16
      t31 = t2 * s26
      t32 = ddreal(5) * s25
      t33 = ddreal(31) * s16
      t34 = ddreal(53) * t11
      t35 = ddreal(22) * t13
      t36 = ddreal(85) * t12
      t37 = s16 + s25
      t38 = ddreal(37) * t19 * t37
      t39 = ddreal(28) * t13
      t40 = t20 * t6 - t38 - t39
      t41 = s25 * (-ddreal(9) * s16 + t16)
      t42 = ddreal(23)
      t43 = ddreal(6)
      t44 = ddreal(21) * s16
      t45 = t43 * t11
      t46 = ddreal(60) * t13
      t47 = ddreal(39)
      t48 = ddreal(10) * s25
      t49 = ((s16 * t47 - t48) * s25 + ddreal(66) * t11) * s25
      t50 = ddreal(80) * t13
      t51 = ddreal(17) * t19
      t52 = -ddreal(12) * t18
      t53 = t51 - t52
      t54 = t2 * s16
      t55 = -t54 - t16
      t56 = ddreal(37) * s16
      t57 = (-t44 + t48) * s25
      t58 = -ddreal(35) * t11
      t59 = ddreal(4)
      t60 = t59 * s25
      t61 = s16 + t60
      t62 = -ddreal(49)
      t63 = ddreal(20) * s25
      t64 = ddreal(15)
      t65 = t59 * s16
      t66 = t64 * s25
      t67 = t11 * t26
      t68 = ddreal(30) * s16
      t69 = ddreal(9) * s25
      t70 = ddreal(16) * s16
      t71 = ddreal(25) * t11
      t72 = ddreal(3) * t28
      t73 = ddreal(3) * s26
      t74 = (t1)**2
      t75 = t11 * t6
      t76 = (s16 - t16) * t74
      t77 = s16 * t4
      t78 = t74 * (s25 * t10 + t15)
      t79 = ddreal(27) * s16
      t80 = t12 * ((-ddreal(21) * s25 + t79) * s25 - ddreal(9) * t11)
      t81 = t2 * s25
      t82 = ddreal(9) * t13
      t83 = ddreal(9) * s26
      t84 = t20 * s25
      t85 = t43 * s16
      t86 = t13 * s25
      t87 = t86 * t1
      t88 = (s34)**2
      t89 = (t88)**2
      t90 = t88 * t89
      t91 = s34 * t89
      t92 = s34 * t88
      t93 = ddreal(7) * s16
      t94 = ddreal(5) * t11
      t95 = s16 * t74
      t38 = ((((t55 + s34) * s34 + t41) * s34 + t25) * s34 - t12 + t19 *
     & ((t16 - t93) * s25 + t94)) * s56 - t14 * t43 - t2 * t91 + ddreal(
     &21) * t87 + ddreal(3) * s16 * (t6 * (s16 + s25) - t76 * s26) + (((
     &(t54 + t73 + t84) * s34 + (-t69 - t85) * s26 + t51 - t52) * s34 - 
     &((-t83 - t84) * s25 + t79 * s26) * s25 - t38 - t39) * s34 + s25 * 
     &(t3 * (s16 * t64 - t81) - t82) + ddreal(22) * t67 + t73 * t25) * s
     &34
      t38 = t38 * s56 + t43 * s16 * (-t78 * s26 + t77) - ddreal(3) * s16
     & * (t76 * t28 + t75) - (((((-t59 * (s16 - s26) + t32 - s34) * s34 
     &- (t59 * t61 + t73) * s26 - t57 - t58) * s34 - s26 * (ddreal(3) * 
     &s26 * t55 + t2 * t53) - t49 - t50) * s34 - t2 * t40 * s26 + t35 * 
     &s25 + t34 * t3 + t33 * t6 + t36 - ddreal(5) * t4 + t69 * (-s25 + t
     &9) * t28) * s34 - s16 * ((((t69 + t70) * s25 + t71) * s25 - ddreal
     &(45) * t13) * s25 + ddreal(44) * t12) - (s25 * (t3 * (-t60 + t68) 
     &- ddreal(18) * t13) + ddreal(44) * t67) * s26 + t5 - t72 * t25) * 
     &s34 + t80
      t39 = ddreal(7) * s25
      t51 = ddreal(11) * s16
      t52 = t59 * t11
      t67 = t88 + t11 + t3
      t96 = -t2 * (s34 * t37 + t19) + t67
      t97 = (-s34 + t37 + s26) * s26 + (s16 + s25 - s34 + t31 + s56) * s
     &56 + t19
      t98 = s25 * t37
      t99 = ddreal(5) * s16
      t100 = (-t81 - t99) * s25 + t52
      t101 = s16 - t81
      t102 = t43 * s25
      t103 = ddreal(10) * t11
      t104 = ddreal(13) * t13
      t105 = t19 * t2 + ddreal(3) * t26
      t106 = ddreal(3) * t11
      t107 = ddreal(5) * t3
      t108 = ((t63 + t33) * s25 + ddreal(18) * t11) * s25
      t109 = ddreal(50) * t13
      t110 = ddreal(22) * s25
      t111 = s16 * t42
      t112 = ddreal(24) * t11
      t113 = (t110 + t111) * s25 + t112
      t114 = s16 + t81
      t115 = t13 + t6
      t116 = ddreal(13) * s25
      t117 = ddreal(12) * s16
      t118 = ddreal(3) * t29
      t119 = ddreal(46) * s16
      t120 = (ddreal(30) * s25 + t119) * s25
      t121 = ddreal(50) * t11
      t122 = ddreal(40) * t13
      t123 = ddreal(25) * s16
      t124 = t59 * s26
      t125 = ddreal(12) * s25
      t126 = ddreal(3) * t105 * s26
      t127 = ddreal(32) * s16
      t128 = t20 * s26
      t129 = s16 * t6
      t130 = t73 * t101
      t131 = t74 * s26
      t132 = ddreal(14) * t12
      t133 = ddreal(10) * s16
      t134 = ddreal(16) * t13
      t135 = t89 + t12
      t110 = ((s25 * ((-t81 + t99) * s25 - t52) - (t105 + t88) * s34 + t
     &13 + ddreal(3) * t114 * t88) * s56 + s25 * ((ddreal(12) * t11 - t7
     &) * s25 - t104) + t59 * t135 - (((ddreal(18) * s25 + t73 + t70) * 
     &s34 - s16 * (t30 + t83) - (ddreal(18) * s26 + t110 + t111) * s25) 
     &* s34 + ((t133 + t102) * s25 - t11 * t20) * s25 + t126 + t134) * s
     &34 - t129 + t130 * t74) * s56 + t59 * t129 * t1 - (((-(t128 + t123
     & + t63) * s34 + (ddreal(36) * s25 + t73 + t127) * s26 + t120 + t12
     &1) * s34 - (t2 * t113 + t83 * t114) * s26 - t108 - t109) * s34 + (
     &((t32 + t85) * s25 + t103) * s25 - t35) * s25 + (t59 * (((t16 + t9
     &9) * s25 - t52) * s25 + t13 * t20) + t126) * s26 + ddreal(25) * t1
     &2) * s34 - ddreal(5) * t91 + ddreal(5) * t14 + t131 * (t100 * t2 +
     & t130) + t82 * t3 - t132 * s25
      t126 = ddreal(3) * t3
      t130 = -t81 - t93
      t136 = s25 * t130 + t94
      t10 = t2 * t11 * t1 - t10 * t3
      t137 = t10 * t28
      t138 = t43 * t3
      t139 = ddreal(11) * s25
      t140 = ddreal(40) * s16
      t141 = ddreal(69) * t12
      t142 = ddreal(16) * t11
      t143 = (t32 + t65) * s25 + t45
      t144 = t143 * s26
      t145 = (((ddreal(58) * s16 + t32) * s25 + ddreal(109) * t11) * s25
     & + ddreal(124) * t13) * s25
      t146 = ddreal(155) * t12
      t147 = ddreal(57) * t11
      t148 = ((t125 + t119) * s25 + t147) * s25 + t46
      t149 = s16 * t37
      t150 = ddreal(38) * t149
      t151 = ddreal(21) * t3
      t152 = t150 + t151
      t153 = t59 * t37
      t154 = ((ddreal(76) * s16 + t48) * s25 + ddreal(131) * t11) * s25
      t155 = ddreal(125) * t13
      t156 = ddreal(75) * t149
      t157 = ddreal(26) * t3
      t158 = t156 + t157
      t159 = ddreal(17) * s16
      t160 = t159 + t116
      t161 = ddreal(14) * s25
      t162 = (ddreal(50) * s16 + t48) * s25
      t163 = ddreal(16) * s25
      t164 = ddreal(13) * s16
      t165 = ddreal(14) * s16
      t128 = t128 * t10
      t166 = (((ddreal(55) * s16 + t139) * s25 + t147) * s25 + ddreal(32
     &) * t13) * s25
      t167 = ddreal(105) * t12
      t168 = ddreal(20) * t11
      t169 = ((t56 + t66) * s25 + t168) * s25
      t170 = ddreal(42) * t13
      t171 = ddreal(11) * t12
      t172 = t124 * t143
      t173 = ddreal(3) * t169 + ddreal(3) * t170
      t174 = t43 * s26
      t175 = t74 * t28
      t176 = ddreal(13) * t11
      t177 = s26 * t37
      t178 = ddreal(5) * t12
      t68 = (((((s26 + t54) * s34 - s16 * (t69 + t165) - (t73 + t32 + t1
     &64) * s26) * s34 + s16 * ((t140 + t163) * s25 + ddreal(42) * t11) 
     &+ s26 * ((t30 + t161) * s26 + t147 + t162) + t118) * s34 - s16 * (
     &((t161 + t119) * s25 + ddreal(65) * t11) * s25 + ddreal(70) * t13)
     & - (((t160 + s26) * s26 + t158) * s26 + t154 + t155) * s26) * s34 
     &+ s16 * ((((t68 + t102) * s25 + ddreal(44) * t11) * s25 + t122) * 
     &s25 + ddreal(70) * t12) + t28 * (t148 * t2 + t153 * t28) + s26 * (
     &t152 * t28 + t145 + t146)) * s34 - t17 * (s16 * ((((t51 + t125) * 
     &s25 + t142) * s25 - ddreal(5) * t13) * s25 + ddreal(42) * t12) + (
     &(((t127 + t66) * s25 + t142) * s25 + ddreal(36) * t13 + t144) * s2
     &6 + (((t139 + t140) * s25 + t71) * s25 + t134) * s25 + t141) * s26
     & + t5)
      t10 = (t2 * (s34 * t10 + t129) + t59 * (t37 * t92 + t86) - t88 * (
     &t143 + t88) - t12 - t94 * t3) * s56 + (-t131 * t101 * t59 + t2 * t
     &4 - t178) * s16 + ddreal(3) * t75 + ddreal(17) * t87 + ddreal(3) *
     & t91 + (((-(t159 + t124 + t116) * s34 + t150 + t151 + ddreal(16) *
     & t177) * s34 - t169 - t170 - t172) * s34 + (((t60 + t165) * s25 + 
     &t176) * s25 - t35) * s25 + t12 * t42 + t128) * s34
      t10 = t10 * s56 + t43 * s16 * (-t175 * t101 + t77) + t80 - ddreal(
     &3) * (t131 * t136 + t75) * s16 - ddreal(3) * t90 - ((((-(t30 + t83
     & + t161) * s34 + (s25 * t47 + ddreal(51) * s16 + t174) * s26 + t15
     &6 + t157) * s34 - s16 * ((ddreal(114) * s16 + ddreal(92) * s25) * 
     &s25 + ddreal(120) * t11) - ddreal(24) * t28 * t37 - ddreal(24) * t
     &6 - t73 * t152) * s34 + (t144 * t43 + t173) * s26 + t166 + t167) *
     & s34 - ((((ddreal(42) * s16 + t125) * s25 + t11 * t47) * s25 - ddr
     &eal(66) * t13) * s25 + t141) * s26 - t6 * ((t81 + t165) * s25 + t1
     &68) - ddreal(48) * t12 * t1 - ddreal(12) * t13 * t3 - ddreal(12) *
     & t137) * s34
      t10 = t10 * s56 + s16 * (-t175 * (t124 * t101 + ddreal(3) * t136) 
     &- t171 * t3 + ddreal(19) * t14 * s25) - ddreal(7) * t12 * t115 + t
     &43 * t11 * (-t131 * t8 + t77) + ((((((-t32 - t174 - t164 + s34) * 
     &s34 + (ddreal(48) * s16 + ddreal(28) * s25 + t83) * s26 + t147 + t
     &162) * s34 - ((ddreal(3) * t160 + t124) * s26 + t158 * t2) * s26 -
     & t154 - t155) * s34 + ((ddreal(16) * s26 * t37 + ddreal(3) * t152)
     & * s26 + t148 * t59) * s26 + t145 + t146) * s34 - ((((s25 + t111) 
     &* s25 + ddreal(51) * t11) * s25 + ddreal(41) * t13) * s25 + t171) 
     &* s25 - ((t172 + t173) * s26 + t2 * (t166 + t167)) * s26 - ddreal(
     &111) * t14) * s34 + t17 * (((((t60 + t70) * s25 + t45) * s25 + t82
     &) * s25 - ddreal(46) * t12) * s25 + ((((ddreal(34) * s16 + t125) *
     & s25 + t11 * t64) * s25 - t109) * s25 + ddreal(53) * t12 + t128) *
     & s26 + ddreal(43) * t14)) * s34
      t10 = -t2 * t12 * (t3 * (-t3 + t11) + t12) - ddreal(3) * t11 * (t1
     &3 * t6 + t175 * t8) + (-t74 * t101 * (t28)**2 - t74 * t136 * t29) 
     &* s16 + s34 * (t68 * s34 + (t17)**2 * (s16 * (((t52 + t138) * s25 
     &- t134) * s25 + t132) + ((((t60 + t133) * s25 + t11) * s25 - ddrea
     &l(14) * t13) * s25 + t12 * t64) * s26 - t2 * (t18 * t6 - t137))) +
     & t10 * s56 - t131 * t13 * ((-t99 - t102) * s25 + t21) + t32 * t13 
     &* t12
      t37 = t20 * t115
      t47 = (s16 + t139) * s25 + t106
      t64 = ddreal(11) * t19
      t26 = ddreal(12) * t26
      t68 = t9 + t39
      t71 = t126 * (s16 - s25 + s34) * s34 * (s56)**3
      t7 = s25 * s34 * (((((-t31 - t99 - t60 + s34) * s34 + t177 * t20 +
     & t138 + ddreal(10) * t149 + t28) * s34 - ((t60 + t93) * s25 + t45)
     & * s25 - ddreal(10) * t13 - s26 * (s26 * t68 + t26 + t64)) * s34 +
     & (t65 * t3 + t118 - t24) * s25 + t178 + (t47 * s26 + t19 * (-t54 +
     & t39) + t37) * s26 + t27) * s34 + (-((-t83 * t1 + t107 - t45) * s2
     &5 - ((-t68 + s34) * s34 + (s16 + t83 + t139) * s25 + t106) * s34 +
     & t13) * s56 - s25 * (s16 * ((-t9 + t60) * s25 - t94) - ddreal(9) *
     & t1 * t28) - t2 * (((-t45 + t107) * s25 + t13) * s26 + t92 * (-t15
     &3 - s26 + s34) + t12 + t4) + s34 * (t2 * t47 * s26 - t15 * s25 - (
     &(t161 + t85) * s26 + t26 + t64) * s34 + t69 * t28 + t93 * t3 + t37
     &)) * s56 - t1 * t17 * ((-t7 - t72) * s25 - s16 * (-t11 + t98) + ((
     &-t54 - t32) * s25 + t11) * s26)) + t71
      t26 = s16 + s25 + s26 - s34 + s56
      t15 = (-s25 - t65) * s25 + t15
      t37 = -s25 + t99
      t47 = (-t16 + t70) * s25 + t45
      t64 = s25 * s26
      t54 = -s25 + t54
      t6 = (t2 * s16 * (-t13 + t6) + ddreal(3) * (s25 * ((ddreal(10) * s
     &26 + t9 + t60) * s16 - t3) - t2 * (s26 * t3 - t13)) * s34 + ddreal
     &(3) * t131 * s25 - t88 * (-t54 * s34 + (t163 + t85) * s16 - ddreal
     &(3) * s25 * (s25 + s26)) + t4 - t11 * s25 * (t69 - t23) + s25 * (-
     &t2 * s25 * (s16 + s34) + t67 + t133 * s34) * s56) * s56 - s25 * (t
     &13 * (-t164 + t125) - t72 * t74) + t2 * (-t131 * t15 + t77) - t59 
     &* s16 * t135 + ((-t2 * t47 * s26 - t176 * s25 + t72 * s25 + ((t124
     & + t116 + t70) * s16 - t31 * s25) * s34 - t117 * t3 - ddreal(24) *
     & t13) * s34 + s16 * (((ddreal(20) * s16 + s25) * s25 - t176) * s25
     & + t134) + t43 * s26 * (((-s25 + t65) * s25 + t106) * s25 + t24 + 
     &t64 * t37)) * s34 + t75
      t9 = s16 + s26 - s34 + s56
      t24 = ddreal(1) / (t96)**2
      t60 = ddreal(1) / (t97)**2
      t4 = s25 * t26 * (-t2 * s16 * (t89 * (-t31 - t16 - t99 + s34) - t1
     &4) - (((s16 * ((t117 + t102) * s25 + t168) + s26 * (s16 * (t116 + 
     &t70) + s26 * t54)) * s34 - s16 * (-t130 * t3 + ddreal(20) * t13) +
     & (-s16 * ((t164 + t125) * s25 + t112) + (-t47 + t64) * s26) * s26)
     & * s34 + t17 * ((((-t16 + t165) * s25 - t11) * s25 + t13 * t43) * 
     &s26 + t22 + t45 * t3 + t129 * t59 - ddreal(12) * t86 + t81 * t37 *
     & t28)) * s34 - t6 * s56 - t11 * t4 - t74 * s25 * t29 + t175 * t15 
     &+ t95 * t100 * s26 + t12 * s25 * (t32 - t85)) * MB1110(1)
      t4 = -t10 * MB1011(1) - t14 * ((-t39 + t51) * s25 - t52) + ddreal(
     &3) * s16 * (t95 * t8 * s26 + t75 * (s16 - s25) + t78 * t28) - ((((
     &((s26 - t65) * s34 + s16 * (t30 + t66) + (-t31 - t32 + t65) * s26)
     & * s34 + s16 * ((s16 * t62 - t63) * s25 - ddreal(60) * t11) + s26 
     &* (t57 + t58 + t28) + t2 * t61 * t28) * s34 + s16 * (((t56 + t48) 
     &* s25 + ddreal(46) * t11) * s25 + t50) + ((s26 * t55 + t53) * s26 
     &+ t49 + t50) * s26) * s34 + t11 * (((-s25 * t42 - t44) * s25 + t45
     &) * s25 - t46) - ((-s26 * t41 - t40) * s26 + (((-t32 + t33) * s25 
     &+ t34) * s25 + t35) * s25 + t36) * s26) * s34 + t17 * (ddreal(11) 
     &* t27 * s16 + t12 * (-ddreal(29) * s25 + t30) + t25 * t28 - t5 + t
     &31 * (-t18 * t19 * t20 + t3 * (-t3 + t21) + t22))) * s34 - t38 * s
     &56 + t76 * s16 * t29 - t4
      result = -t24 * t60 * ((s16 * I300s16(1) - MB1001(0)) * (t2 * (t11
     & * (t18 * t3 + t12) + t90) + ((((-(ddreal(5) * s26 + t84 + t117) *
     & s34 + (t79 + t125) * s25 + (t123 + t124 + t63) * s26 + ddreal(30)
     & * t11) * s34 - ((t84 + t111) * s25 + ddreal(28) * t11) * s25 - s2
     &6 * (t120 + t121 + t28) - t122 - t2 * (t69 + t23) * t28) * s34 + t
     &2 * s25 * t115 + s16 * (t3 * (t116 + t117) + ddreal(30) * t13) + s
     &26 * (s26 * t113 + t108 + t109) + t118 * t114) * s34 - t17 * (ddre
     &al(12) * t13 * t1 + t3 * (t106 + t107) + s26 * (t105 * s26 + ((t10
     &2 + t93) * s25 - t103) * s25 + t104))) * s34 + t110 * s56 + t131 *
     & (s16 * (-t59 * t98 + t94) + (s26 * t101 + t100) * s26) + t86 * (t
     &126 - t94)) + t4 / s16) / ddreal(4) - t60 * (-t7 * MB1101(0) * t24
     & + (t9)**2 * t26 * MB1111D2(0)) / ddreal(2)

           intHs16s25s26s34s56x1213D6eps0 = result
       end function intHs16s25s26s34s56x1213D6eps0

       function intHs16s25s26s34s56x1213D6eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1213D6eps1

           type(dd_complex) :: result

      result = ddreal(0)

           intHs16s25s26s34s56x1213D6eps1 = result
       end function intHs16s25s26s34s56x1213D6eps1

       function intHs16s25s26s34s56x1220D4eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1220D4eps0
           type(dd_complex) ::  t1,t2

           type(dd_complex) :: result

      t1 = s16 + s25 + s26 - s34 + s56
      t2 = s16 + s26 - s34 + s56
      t1 = ddreal(1) / (t1)**2
      result = -t1 * (((s16 + s26 + s56) * MB1011(1) - MB1001(0) - ddrea
     &l(1)) * (s16 - s25 + s26 - s34 + s56) + s25 * (s25 - ddreal(5) * t
     &2) * MB1110(1)) / ddreal(2) + s25 * t2 * MB1110(0) * t1

           intHs16s25s26s34s56x1220D4eps0 = result
       end function intHs16s25s26s34s56x1220D4eps0

       function intHs16s25s26s34s56x1220D4eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1220D4eps1
           type(dd_complex) ::  t1,t2

           type(dd_complex) :: result

      t1 = s16 + s25 + s26 - s34 + s56
      t2 = ddreal(1) / ddreal(2)
      t1 = ddreal(1) / (t1)**2
      result = t1 * (t2 * (s16 - s25 + s26 - s34 + s56) + s25 * (s16 + s
     &26 - s34 + s56) * MB1110(1))

           intHs16s25s26s34s56x1220D4eps1 = result
       end function intHs16s25s26s34s56x1220D4eps1

       function intHs16s25s26s34s56x1221D4eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1221D4eps0
           type(dd_complex) ::  t1,t10,t11,t12,t13,t2,t3,t4,t5,t6,t7,t8
           type(dd_complex) ::  t9

           type(dd_complex) :: result

      t1 = s16 + s26 + s56
      t2 = ddreal(2) * s26
      t3 = s16 * s25
      t4 = (-s34 + s16 + s25 + s26) * s26 + (t2 + s16 + s25 - s34 + s56)
     & * s56 + t3
      t5 = s16 + s26 - s34 + s56
      t6 = s16 + s25 + s26 - s34 + s56
      t7 = s16 + s26
      t8 = (s26)**2
      t9 = s16 - s25
      t10 = t9 * s26
      t4 = ddreal(1) / t4
      t6 = ddreal(1) / t6
      t11 = ddreal(1) / s16
      t12 = ddreal(1) / t5
      t13 = t4 * t11
      result = ddreal(2) * t13 * t1 - t6 * (t11 * (-s25 * MB1110(0) + (-
     &s25 * (s25 * s26 + (s16 - s26) * s34 + (t2 + s25 - s34 + s56) * s5
     &6 + t3 - (s16)**2 + t8) * MB1110(1) - t1 * (t7 * s16 + (-ddreal(2)
     & * s16 - s25 - s26 + s34) * s34 + (s16 - s34) * s56) * MB1011(1)) 
     &* t4) + t5 * MB1001(0) * t4) + t1 * I300s16(1) * t4 + t13 * (((s16
     & * t9 + (s34)**2) * s16 + t9 * ((s56)**2 + t8) + ddreal(2) * (-s34
     & * t7 + t10) * s16 + ddreal(2) * ((-s25 - s34 + s16) * s16 + t10) 
     &* s56) * MB1111D2(0) + s25 * s34 * t1 * MB1101(0)) * t12

           intHs16s25s26s34s56x1221D4eps0 = result
       end function intHs16s25s26s34s56x1221D4eps0

       function intHs16s25s26s34s56x1221D4eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1221D4eps1
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s16 + s25 + s26 - s34 + s56
      t1 = ddreal(1) / t1
      result = ddreal(1) / s16 * t1 * (s25 * MB1110(1) + ddreal(1))

           intHs16s25s26s34s56x1221D4eps1 = result
       end function intHs16s25s26s34s56x1221D4eps1

       function intHs16s25s26s34s56x1222D6eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1222D6eps0
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           type(dd_complex) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           type(dd_complex) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           type(dd_complex) ::  t42,t43,t44,t45,t46,t47,t48,t49,t5,t50,t51,t52
           type(dd_complex) ::  t53,t54,t55,t56,t57,t58,t59,t6,t60,t61,t62,t63
           type(dd_complex) ::  t64,t65,t66,t67,t68,t69,t7,t70,t71,t72,t73,t8
           type(dd_complex) ::  t9

           type(dd_complex) :: result

      t1 = ddreal(2)
      t2 = t1 * s25
      t3 = ddreal(7) * s16
      t4 = s16 - s25
      t5 = t1 * s16
      t6 = -s25 + t5
      t7 = s16 + s26
      t8 = (s16)**2
      t9 = (t8)**2
      t10 = s16 * t8
      t11 = s25 * t4
      t12 = s16 + s25
      t13 = (s25)**2
      t14 = s25 * t13
      t15 = s16 * t13
      t16 = ddreal(11)
      t17 = ddreal(4)
      t18 = (s26)**2
      t19 = ddreal(3)
      t20 = ddreal(18)
      t21 = ddreal(6) * s16
      t22 = t19 * s25
      t23 = ddreal(5) * s25
      t24 = ddreal(5) * s16
      t25 = (t24 + t2) * s25
      t26 = t1 * s26
      t27 = t17 * s16
      t28 = ddreal(7) * s25
      t29 = ddreal(16) * s16
      t30 = s16 + t2
      t31 = t12 * s26
      t32 = ddreal(7) * t8
      t33 = t18 + t13
      t34 = ddreal(6) * s25
      t35 = t17 * s26
      t36 = t19 * t33
      t37 = (t4)**2
      t38 = s16 * t14
      t39 = t37 * s26
      t40 = s25 * t16
      t41 = t37 * t18
      t42 = t41 * t19
      t43 = t19 * t31
      t44 = t19 * s26
      t45 = (s34)**2
      t46 = (t45)**2
      t47 = s34 * t45
      t48 = t17 * s25
      t49 = t39 * t19
      t50 = s16 * s25
      t51 = -t1 * (s34 * t12 + t50) + t13 + t45 + t8
      t37 = t37 * s26 * t18
      t52 = t8 * t14
      t53 = t39 * s16
      t54 = t17 * t8
      t6 = t10 * ((-ddreal(9) * s16 + t34) * s25 + t54) + ((((-s26 + t27
     &) * s34 - s16 * (t29 + t28) + (t22 - t27 + t26) * s26) * s34 + s16
     & * (ddreal(24) * t8 + t25) - ((t23 + t5 + s26) * s26 - (-t22 + t21
     &) * s25 - t20 * t8) * s26) * s34 - t7 * ((-t16 * t8 - t13) * s25 +
     & ddreal(16) * t10 + t17 * (-(-t11 - t8) * s26 + t15) - t1 * t12 * 
     &t18)) * s34 + ((-t51 * s56 + t1 * ((-s16 * t12 + t1 * t13 + t43 + 
     &t45) * s34 + t10) - (t23 + t44 + t5) * t45 - t14 - t49 + t50 * (-t
     &24 + t48)) * s56 - t1 * (-t39 * t6 + t38) + t8 * ((-t29 + t40) * s
     &25 + t32) - (((t17 * (s16 - s26) - t22 + s34) * s34 - (s16 * t20 +
     & t34 - t35) * s16 + t36 + ddreal(10) * s25 * s26) * s34 - (t32 + t
     &13) * s25 - (-t4 * t30 * t17 + ddreal(6) * t31) * s26 + ddreal(20)
     & * t10) * s34 - t42) * s56 - t37 - t52 + t53 * (t3 - t2) + t41 * t
     &6
      t20 = (-s34 + t12 + s26) * s26 + (t26 + s16 + s25 - s34 + s56) * s
     &56 + t50
      t29 = s16 * t7
      t55 = ddreal(6) * t29
      t33 = (t44 + t5) * s25 + t33 + t55
      t56 = s25 + s26
      t57 = t8 + t13
      t58 = s16 * t57
      t59 = t54 * s25
      t11 = t1 * (t58 + t39 - t47) - s34 * (-(t22 + t26 + t21) * s34 + (
     &t27 + t2) * s26 - t11 + ddreal(6) * t8) - t59
      t60 = t7 * ((t27 - t2) * s16 + (s25 + t5) * s26 + t13)
      t61 = (-t1 * s16 * (s25 + s34) + (s34 - s25) * s34 + t13 + t8) * s
     &56
      t62 = s25 * t8 - t39
      t63 = t62 * s16 + t47 * (t5 + t56)
      t64 = t57 * t8
      t65 = s25 + t24
      t66 = ddreal(10) * t8
      t67 = s25 + t27
      t68 = t67 * t18
      t69 = t19 * s16
      t70 = ddreal(20) * s16
      t71 = ddreal(13) * s16
      t72 = ddreal(15) * s26
      t73 = ddreal(8) * s26
      t16 = (-(t1 * t58 - t59 - (s16 * t65 + (-t67 + s34) * s34) * s34) 
     &* s56 - t1 * t46 - ddreal(6) * (t58 + t39) * s16 + (s16 * (t19 * s
     &25 * t56 + (-t23 + t72 + t70) * s16) + ((ddreal(12) * s16 + t44 + 
     &t48) * s34 - s16 * (ddreal(24) * s16 + ddreal(12) * s26) - (s16 * 
     &t16 + t2 + t44) * s25) * s34) * s34 + ddreal(12) * t10 * s25) * s5
     &6 + ddreal(12) * t8 * t62 - ddreal(6) * (t64 + t41) * s16 + ((((-d
     &dreal(10) * s16 - t22 - t35 + s34) * s34 + (ddreal(30) * s16 + ddr
     &eal(17) * s25 + ddreal(24) * s26) * s16 + t36 + t73 * s25) * s34 -
     & s16 * ((t28 + t71) * s25 + ddreal(40) * t8) - ((ddreal(22) * s16 
     &+ t48) * s25 + ddreal(48) * t8) * s26 - t14 - t68 * t19) * s34 + t
     &29 * (s16 * (ddreal(25) * s16 + t72) + (t44 + t34 - t71) * s25)) *
     & s34
      t34 = t8 + t45
      t36 = -s16 + s25 + s34
      t47 = t10 - t47
      t56 = (-s25 + t69 + s26) * s26
      t58 = ddreal(8) * t10
      t59 = ddreal(15) * s16
      t62 = s25 * t12
      t71 = s16 * s34
      t72 = s25 * s34
      t12 = t1 * t50 * (t9 + t46) + s25 * (((-t71 * (ddreal(8) * s16 + d
     &dreal(5) * s26 + t48) + s16 * (t25 + ddreal(12) * t8) + s26 * (s16
     & * (t28 + t59) + s26 * t67)) * s34 - t7 * (-t1 * t50 * t4 + s26 * 
     &(s25 * t30 + t31 + t32) + t58)) * s34 + (((-t50 * t1 + t13 - t71 -
     & t72 + t8) * s56 + s16 * ((-t3 + t2) * s25 + t54) + t49 + (s34 * t
     &67 - t1 * t62 - t43 - ddreal(8) * t8) * s34 + t14) * s56 + ddreal(
     &5) * s16 * t47 + t1 * (t39 * t67 + t38) + t42 + (((t73 + t28 + t59
     &) * s16 + t26 * s25) * s34 - t17 * ((t54 + t62) * s26 + t15) - t8 
     &* (-s25 + t59) - t19 * t12 * t18) * s34 + t13 * t8 - t58 * s25) * 
     &s56 + t37 + t52 + t53 * (t24 + t2) + t41 * t67) - t19 * t9 * t13
      t25 = ddreal(1) / t51
      t20 = ddreal(1) / (t20)**2
      t27 = (s16 + s26 - s34 + s56) * (t1 * t34 - (t44 + t27 + t2) * s34
     & + (t19 * (s16 - s34) - s25 + t26 + s56) * s56 - t50 + t56) * MB11
     &11D2(0)
      t11 = t25 * t20 * ((-t1 * t63 * s16 + s16 * (((t45 + t33) * s34 - 
     &t60) * s34 + (t61 + t11) * s56 + t8 * (t8 + t13) + t41)) * I300s16
     &(1) + (t1 * t63 - ((t45 + t33) * s34 - t60) * s34 - (t61 + t11) * 
     &s56 - t41 - t64) * MB1001(0))
      result = t20 * (t25 * (((t1 * (t10 * t57 + t37) * s16 - s34 * ((((
     &(s26 + t5) * s34 - (t22 + t26) * s26 - ddreal(10) * t29 - ddreal(5
     &) * t50) * s34 + t17 * (t18 * (s25 + t69) + t15) + t8 * (t40 + t70
     &) + s26 * (t18 + (ddreal(17) * s16 + t22) * s25 + ddreal(30) * t8)
     &) * s34 - t7 * (s16 * ((t23 + t69) * s25 + ddreal(20) * t8) + t14 
     &+ t68 + t26 * ((s25 + t24) * s25 + t66))) * s34 + s16 * (t7)**2 * 
     &((-t3 + t22) * s25 + s26 * t65 + t66)) - t16 * s56 + ddreal(6) * t
     &39 * t8 * t7 - t48 * s16 * t9) * MB1011(1) - t12 * MB1110(1) + t6)
     & / s16 - t72 * (-t1 * t47 - ((t44 + t48 + t21) * s34 - (t35 + t2) 
     &* s25 - t18 - t50 - t55) * s34 - (-t36 * s56 - t1 * t36 * s26 - t1
     &7 * s25 * (s16 + s34) + t19 * t34 + t13 - t21 * s34) * s56 - t15 -
     & t56 * t4 + t22 * t8) * MB1101(0)) + t27) / ddreal(2) + t11

           intHs16s25s26s34s56x1222D6eps0 = result
       end function intHs16s25s26s34s56x1222D6eps0

       function intHs16s25s26s34s56x1222D6eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1222D6eps1

           type(dd_complex) :: result

      result = ddreal(0)

           intHs16s25s26s34s56x1222D6eps1 = result
       end function intHs16s25s26s34s56x1222D6eps1

       function intHs16s25s26s34s56x1231D6eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1231D6eps0
           type(dd_complex) ::  t1,t10,t100,t101,t102,t103,t11,t12,t13,t14,t15,t16
           type(dd_complex) ::  t17,t18,t19,t2,t20,t21,t22,t23,t24,t25,t26,t27
           type(dd_complex) ::  t28,t29,t3,t30,t31,t32,t33,t34,t35,t36,t37,t38
           type(dd_complex) ::  t39,t4,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49
           type(dd_complex) ::  t5,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t6
           type(dd_complex) ::  t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t7,t70
           type(dd_complex) ::  t71,t72,t73,t74,t75,t76,t77,t78,t79,t8,t80,t81
           type(dd_complex) ::  t82,t83,t84,t85,t86,t87,t88,t89,t9,t90,t91,t92
           type(dd_complex) ::  t93,t94,t95,t96,t97,t98,t99

           type(dd_complex) :: result

      t1 = ddreal(9)
      t2 = (s16)**2
      t3 = (t2)**2
      t4 = s16 * t3
      t5 = s16 * t2
      t6 = ddreal(16) * s25
      t7 = ddreal(37) * s16
      t8 = t2 * t1
      t9 = ddreal(13) * t5
      t10 = ddreal(12)
      t11 = ddreal(8) * s25
      t12 = t10 * t5
      t13 = ((ddreal(47) * s16 + t11) * s25 + ddreal(35) * t2) * s25
      t14 = t13 - t12
      t15 = ddreal(2)
      t16 = t15 * s16
      t17 = ddreal(19)
      t18 = t17 * s25
      t19 = ddreal(8) * s16
      t20 = ddreal(14) * s25
      t21 = t20 + t19
      t22 = s16 + s26
      t23 = (s25)**2
      t24 = s25 * t23
      t25 = -t2 + t23
      t26 = -ddreal(22)
      t27 = ddreal(18) * s16
      t28 = ddreal(4)
      t29 = t28 * s16
      t30 = -ddreal(31) * s25
      t31 = (s26)**2
      t32 = (t31)**2
      t33 = s26 * t31
      t34 = s16 * t1
      t35 = ddreal(7) * s25
      t36 = ddreal(16) * t2
      t37 = ddreal(3)
      t38 = t37 * s16
      t39 = t38 - s25
      t40 = t29 + s25
      t41 = t16 * t23
      t42 = t37 * s26
      t43 = t2 * s25
      t44 = t37 * s25
      t45 = ddreal(6) * s26
      t46 = ddreal(10) * s16
      t47 = ddreal(16) * s16
      t48 = -t29 + s26
      t49 = ddreal(42) * s16
      t50 = -ddreal(93) * s25
      t51 = -ddreal(29)
      t52 = ddreal(38) * t2
      t53 = s16 + s25
      t54 = s25 * t53
      t55 = ddreal(40) * s26
      t56 = t10 * s26
      t57 = ddreal(15) * s26
      t58 = (t16 + s25) * (t18 + s16)
      t59 = t58 * t37
      t60 = ddreal(60)
      t61 = ddreal(5) * s25
      t62 = t1 * s26
      t63 = (s34)**2
      t64 = (t63)**2
      t65 = s34 * t63
      t66 = ddreal(14) * s16
      t67 = ddreal(32) * s16
      t13 = (((-ddreal(10) * s34 + t19 + t20 + t57) * s56 + (ddreal(39) 
     &* s16 + t18) * s25 + (ddreal(56) * s25 + ddreal(30) * s26 + t67) *
     & s26 + (s34 * t10 + t30 - t55 - t66) * s34 + t15 * t2) * s56 + (dd
     &real(30) * t31 + t59) * s26 + (t28 * (t61 + t62) - ddreal(6) * s34
     &) * t63 - t12 + t13 + ((s25 * t26 - ddreal(49) * s16) * s25 + (-s2
     &6 * t60 - t49 + t50) * s26 + ddreal(18) * t2) * s34 + ddreal(6) * 
     &t31 * t21) * s56 + s16 * (((t7 + t6) * s25 + t8) * s25 - t9) + t31
     & * ((ddreal(4) * t21 + t57) * s26 + t59) + t15 * t14 * s26 + ((t55
     & * s25 + (-t44 + t46 - t56 + s34) * s34 + t37 * t54 - ddreal(36) *
     & t2 + ddreal(36) * t31) * s34 + s16 * ((s25 * t51 - t34) * s25 + t
     &52) - ddreal(40) * t33 + ((t50 - t49) * s26 + (-ddreal(98) * s16 -
     & ddreal(44) * s25) * s25 + ddreal(36) * t2) * s26 - t24) * s34
      t18 = t3 * s25
      t20 = (s56)**2
      t49 = (t20)**2
      t50 = s56 * t49
      t51 = s56 * t20
      t55 = s16 + s25 + s26 - s34 + s56
      t57 = t15 * s26
      t59 = s16 * s25
      t68 = (-s34 + t53 + s26) * s26 + (t57 + s16 + s25 - s34 + s56) * s
     &56 + t59
      t69 = t15 * s25
      t70 = -s16 + s34
      t71 = t70 * s16
      t72 = s25 + s26
      t73 = s34 * t22
      t74 = t2 + t63
      t75 = s16 + s26 - s34 + s56
      t76 = (-t44 - s16) * s25
      t77 = ddreal(10) * s26
      t78 = ddreal(6) * s16
      t79 = ddreal(5) * s26
      t80 = t28 * s25
      t81 = t28 * s26
      t82 = s34 * t15
      t83 = t5 - t65
      t84 = (-s25 * t1 - t38) * s25 + t36
      t85 = ddreal(5) * s16
      t86 = (-t44 - t85) * s25 + ddreal(14) * t2
      t87 = ddreal(7) * s16
      t88 = -t69 + t87
      t89 = ddreal(8) * t2
      t90 = ddreal(20) * t2
      t91 = ddreal(20) * s16
      t92 = ddreal(42) * t2
      t93 = (t91 + t44) * s25
      t94 = ddreal(25) * s16
      t95 = t16 + s26
      t96 = ddreal(6) * s25
      t97 = ddreal(13) * s16
      t98 = ddreal(26)
      t60 = t60 * s16
      t99 = ddreal(8) * t39
      t100 = s16 * s26
      t101 = t1 * t23
      t102 = s16 * t23
      t103 = (t22)**2
      t54 = (-(t15 * s25 * t70 + (-ddreal(7) * s34 + t79 + t78) * s16 + 
     &t63) * s56 - s16 * ((t66 - t61) * s16 + (t99 + t77) * s26) + t37 *
     & (t102 + t65) - s34 * ((t91 + t81 - s25) * s34 + (-ddreal(28) * s1
     &6 + t11) * s26 + t28 * t54 - ddreal(31) * t2)) * s56 - s16 * ((-t1
     &01 + t36) * s16 + t31 * (t10 * t39 + t77)) - t37 * ((s26 * t86 - t
     &43) * s16 + t64) - ((-(t94 + t80 + t62) * s34 + (-s25 + t97) * s25
     & + (-t44 + t45 + t60) * s26 + ddreal(57) * t2) * s34 + ((ddreal(11
     &) * s16 + t69) * s25 - ddreal(6) * t2) * s25 - ddreal(51) * t5 + s
     &26 * (-ddreal(6) * s26 * t88 + t10 * t54 - ddreal(93) * t2)) * s34
      t91 = t103 * (((-t69 - t38) * s25 + t89) * s25 + ddreal(10) * t5 +
     & s26 * (s26 * t88 - t23 * t28 + ddreal(17) * t2)) * s34
      t17 = t22 * (s16 * ((t80 + t27) * s25 + t90) + ((s16 * t17 - s25 +
     & s26) * s26 + (t66 - s25) * s25 + t52) * s26 + t24) * t63
      t52 = s16 * (t32 * t39 + t4)
      t66 = t37 * t23
      t1 = t15 * t52 - t95 * s34 * t64 - (t33 * (-t31 - t86) - t50) * s1
     &6 + s25 * t4 - (t54 * s56 + t1 * t5 * t25 - ((((t44 + t45 + t97 - 
     &s34) * s34 - (ddreal(50) * s16 + t11 + t62) * s26 - t92 - t93) * s
     &34 + s16 * ((t44 + t67) * s25 + ddreal(58) * t2) + t28 * t33 + ((-
     &t44 + t60) * s26 + (s16 * t98 - t69) * s25 + ddreal(114) * t2) * s
     &26 + t24) * s34 - t22 * (((-t80 - t46) * s25 + t36) * s25 + ddreal
     &(37) * t5 + s26 * (t28 * t88 * s26 + (-t96 + t97) * (t69 + t85))))
     & * s34 - t18 - t100 * (((t79 + t99) * s26 + t37 * t86) * s26 + t16
     & * t84)) * s56 + t64 * (t37 * s26 * t72 + (ddreal(13) * s26 + t61 
     &+ t46) * s16) - t65 * (s16 * ((t47 + t80) * s25 + t90) + t33 * t37
     & + s26 * ((t94 + t80) * s26 + t92 + t93)) - t91 + t17 + t5 * (-t1 
     &* t25 + t59) * s26 + t2 * t84 * t31 - t66 * t3
      t4 = s16 - s25
      t17 = t4 * s26
      t52 = t2 + t23
      t54 = s25 * t52
      t59 = (t2 * t28 + t66) * s16 + t54
      t60 = (t96 + t85) * s25 + t89
      t66 = t16 + t44
      t19 = t96 + t19
      t67 = t10 * t2
      t84 = s16 * t22
      t64 = (t5 * t53 + t64) * s16
      t65 = t29 * (t15 * t22 + s25) * t65
      t86 = (t4)**2
      t88 = s26 * t86
      t90 = ddreal(1) / s16
      t55 = ddreal(1) / (t55)**2
      t68 = ddreal(1) / (t68)**2
      t12 = s25 * (t15 * t64 + t37 * (-t2 * t23 * t4 + t31 * t59) + s16 
     &* ((t2 * t37 + ddreal(6) * t23) * s25 + ddreal(8) * t5) * s26 + ((
     &t12 + t41 + t37 * t40 * t31 + ddreal(10) * t43 + t100 * (ddreal(24
     &) * s16 + t35)) * s34 - t22 * (t19 * t31 + ddreal(8) * t2 * t53 - 
     &t102 + t57 * (-t76 + t89))) * s34 + (ddreal(6) * s16 * t24 + ddrea
     &l(6) * s26 * t59 + t37 * (s25 * t5 + t31 * t60) + t33 * t28 * t66 
     &+ (((t45 + t87) * s25 + ddreal(24) * t84) * s34 - s16 * ((t61 + t4
     &6) * s25 + ddreal(24) * t2) - (t28 * ((t44 + t29) * s25 + t67) + t
     &45 * t19 / ddreal(2)) * s26) * s34 + ddreal(8) * s16 * t83) * s56 
     &+ t20 * ((t101 + t67) * s16 + t37 * (s26 * t60 + t40 * t63 + t54) 
     &- t82 * ((t44 + t29 + t62) * s25 + t10 * t84) + ddreal(6) * t66 * 
     &t31) + t32 * t66 + t33 * t60 + t49 * t66 + t51 * ((t85 + t56) * s2
     &5 + ddreal(8) * (s26 - s34 + s16) * s16 - ddreal(6) * s25 * (-s25 
     &+ s34)) - t65)
      t19 = t75 * (t15 * t83 + ((-t69 + t29 + s26) * s26 + ddreal(5) * t
     &2 + t76) * s26 + ((t79 + t80 + t78) * s34 - t15 * s25 * t72 - (t78
     & + t61 + t77) * s16 - t28 * t31) * s34 + ((-t28 * t70 + s56 + t42 
     &- t69) * s56 + (-t81 - t82 - s16) * s25 - ddreal(8) * s26 * t70 - 
     &t37 * (t23 - t31) + ddreal(5) * t74 - t46 * s34) * s56 + t43 - t38
     & * t23) * MB1001(0)
      t1 = t68 * ((-t15 * s16 * t74 - t37 * ((s16 * (-t69 + s16) - s25 *
     & s26) * s26 - t43) + ((t37 * t72 - s34 + s56) * s56 + ddreal(6) * 
     &s25 * t22 + t37 * (t71 + t31) - t57 * s34) * s56 + t33 - t73 * t48
     &) * I300s16(1) + t55 * (t90 * (-t1 * MB1011(1) + t2 * (t23 * (t34 
     &+ t11) - t28 * t5) + t37 * (s26 * t32 + t50) + (((s26 * t21 + t58)
     & * s26 + t14) * s26 + s16 * (((t7 + t6) * s25 + t8) * s25 - t9)) *
     & s26 + ((-t42 * t39 * t40 + ddreal(20) * s25 * t31 + (t48 * s34 + 
     &s16 * (t47 + t35) + (-t44 - t45 + t46) * s26) * s34 + t10 * t33 - 
     &t41 - ddreal(15) * t43 - ddreal(24) * t5) * s34 + t22 * (s16 * ((t
     &34 - t35) * s25 + t36) - ddreal(10) * t33 + (-t27 * s25 + (-t29 + 
     &t30) * s26 + t25 * t26) * s26 - t24)) * s34 + t13 * s56 - t18 + t1
     &2 * MB1110(1)) + t19))
      t5 = ddreal(1) / t75
      result = t1 / ddreal(4) + t90 * (((-t37 * (s16 * (t22 * t95 * t63 
     &+ t88 * t22) + ((s16 * ((-t16 + s25 - s26 + s34) * s34 + t2 + t23)
     & + t88 - t43 * t15) * s56 + t15 * s16 * (-t43 + t88) + s16 * (s16 
     &* t52 + ((t57 + t38 - s34) * s34 - t22 * (-t69 + t38 + s26)) * s34
     &) + t31 * t86) * s56) - s16 * ((-t103 * (-t44 + t29 + s26) + t63 *
     & (-t42 - t29 - s25 + s34)) * s34 + t2 * t52) - t33 * t86 - t51 * (
     &-t16 * s25 + t23 - t71) + t69 * t3) * MB1111D2(0) - s25 * s34 * (t
     &15 * ((-t73 + t17) * s16 + ((-s25 - s34 + s16) * s16 + t17) * s56)
     & + (s16 * t4 + t63) * s16 + t4 * (t31 + t20)) * MB1101(0)) * t68 *
     & t5 + t23 * MB1110(0) * t55) / ddreal(2)

           intHs16s25s26s34s56x1231D6eps0 = result
       end function intHs16s25s26s34s56x1231D6eps0

       function intHs16s25s26s34s56x1231D6eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1231D6eps1
           type(dd_complex) ::  t1,t2,t3

           type(dd_complex) :: result

      t1 = ddreal(3)
      t2 = s16 + s25 + s26 - s34 + s56
      t3 = ddreal(1) / ddreal(4)
      t2 = ddreal(1) / (t2)**2
      result = ddreal(1) / s16 * t2 * (t3 * (s25 * t1 + s16 + s26 - s34 
     &+ s56) + (s25)**2 * MB1110(1) / ddreal(2))

           intHs16s25s26s34s56x1231D6eps1 = result
       end function intHs16s25s26s34s56x1231D6eps1

       function intHs16s25s26s34s56x1310D4eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1310D4eps0
           type(dd_complex) ::  t1,t2,t3,t4,t5

           type(dd_complex) :: result

      t1 = s16 + s26 - s34 + s56
      t2 = s16 + s25 + s26 - s34 + s56
      t3 = ddreal(3)
      t4 = ddreal(6)
      t5 = ddreal(4) * s25
      t2 = ddreal(1) / (t2)**2
      result = t2 * (((s16 + s26 + s56) * MB1011(1) - MB1001(0)) * (t1 *
     & t3 + s25) + (t3 * ((s16)**2 + (s26)**2 + (s34)**2 + (s56)**2) + (
     &s16 * t4 - t5) * s26 + (-t4 * (s16 + s26) + t5) * s34 + (t4 * (s16
     & + s26 - s34) - t5) * s56 - (s25)**2 - t5 * s16) * MB1110(1)) / dd
     &real(4) - t1 * t2 * (-t1 * MB1110(0) + ddreal(1)) / ddreal(2)

           intHs16s25s26s34s56x1310D4eps0 = result
       end function intHs16s25s26s34s56x1310D4eps0

       function intHs16s25s26s34s56x1310D4eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1310D4eps1
           type(dd_complex) ::  t1,t2

           type(dd_complex) :: result

      t1 = s16 + s26 - s34 + s56
      t2 = s16 + s25 + s26 - s34 + s56
      t2 = ddreal(1) / (t2)**2
      result = t2 * (-s25 / ddreal(4) - ddreal(3) / ddreal(4) * t1 + (t1
     &)**2 * MB1110(1) / ddreal(2))

           intHs16s25s26s34s56x1310D4eps1 = result
       end function intHs16s25s26s34s56x1310D4eps1

       function intHs16s25s26s34s56x1311D4eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1311D4eps0
           type(dd_complex) ::  t1,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t2
           type(dd_complex) ::  t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t3,t30
           type(dd_complex) ::  t31,t32,t33,t34,t35,t36,t37,t38,t39,t4,t40,t41
           type(dd_complex) ::  t42,t5,t6,t7,t8,t9

           type(dd_complex) :: result

      t1 = s16 - s25
      t2 = s16 + s25
      t3 = ddreal(2)
      t4 = t3 * s25
      t5 = t4 * s26
      t6 = (s34)**2
      t7 = (s25)**2
      t8 = s25 * t7
      t9 = (s16)**2
      t10 = s16 * t9
      t11 = t9 + t6 + t7
      t12 = t9 + t7
      t13 = (t1)**2
      t14 = t13 * s26
      t15 = t12 * s16
      t16 = t4 * t9
      t17 = -t16 + (-s16 * t2 + s26 * s34 - t5) * s34 + (-t3 * s25 * (s1
     &6 + s34) + t11) * s56 + t14 + t15
      t18 = s16 * s25
      t11 = -t3 * (s34 * t2 + t18) + t11
      t19 = s16 + s25 - s34
      t20 = t3 * s26
      t21 = (-s34 + t2 + s26) * s26 + (t20 + t19 + s56) * s56 + t18
      t22 = (s26)**2
      t23 = s16 * t1
      t24 = s26 * t2
      t25 = ddreal(4)
      t26 = t25 * s26
      t27 = (s56)**2
      t28 = t9 * s25
      t29 = s16 + s25 + s26 - s34 + s56
      t30 = s16 + s26 + s56
      t31 = t9 - t7
      t32 = t4 * s16
      t33 = s16 * t7
      t34 = ddreal(3)
      t35 = t34 * s16
      t36 = s16 + t4
      t37 = t3 * s16
      t38 = s25 + t37
      t39 = t23 * t3 + t7
      t40 = t2 * t22
      t41 = s16 * s34
      t21 = ddreal(1) / t21
      t29 = ddreal(1) / t29
      t11 = ddreal(1) / t11
      t42 = ddreal(1) / s16
      t5 = t30 * (-t4 * t10 + (t15 + t14) * s16 + ((-(s16 - s26) * s34 +
     & t34 * t9 - t5 - s26 * (s16 + s26)) * s34 + t19 * t27 + t9 * (-t35
     & + t4) + (t24 + t32 - t31) * s26 + t33) * s34 + (-t16 + ((-t3 * (s
     &25 + s26) - s16 + s34) * s34 + t3 * (t24 + t18) - t31) * s34 + t15
     &) * s56) * MB1011(1)
      t19 = s34 * ((s26 * t1 + t9) * s25 - s26 * t6 + (s26 * t36 + t18) 
     &* s34 - (-s25 * t1 + (-t36 + s34) * s34) * s56 - t33) * MB1101(0)
      t15 = (-t3 * (-s34 * t27 + t28) + ((s16 - t20) * s34 - t3 * (t23 -
     & t22) + t24) * s34 + (-t3 * (t18 + t6) + (t26 + t2) * s34 + t7 + t
     &9) * s56 + t14 + t15) * MB1001(0) * t29
      t24 = (s16 + s26 - s34 + s56) * MB1110(0) * t42 * t29
      t27 = -ddreal(1) / ddreal(2)
      t31 = t17 * t42 * t11 * t21
      result = t27 * (t21 * (t11 * (t17 * I300s16(1) + (t29 * (-(t3 * t2
     &8 * t12 + (-t25 * t28 * t1 - s26 * (t23 * (s25 + t35) + t40) - t3 
     &* t39 * t22) * s34 - ((-(-s25 * s34 - t32 - t41 + t7 + t9) * s56 -
     & t10 * t3 + t34 * (t28 - t14) - (-t34 * s26 * t2 + s34 * t38 - t23
     & * t25 - t3 * t7) * s34 - t8) * s56 - t34 * (s16 * t8 + t13 * t22)
     & - ((-t41 + (s16 * t25 + t4) * s26 + t23 * t34) * s34 - t34 * (t40
     & + t10) - t26 * t39 + t18 * (s25 + t37)) * s34 - t10 * t2 + ddreal
     &(5) * t9 * t7 - t20 * t13 * t38) * s56 + t6 * ((s26 * t38 + t23 * 
     &t34) * s26 + t16) - s16 * s26 * s34 * t6 - t14 * (-s16 * (s25 * t3
     &4 + s16) + (-t38 - s26) * s26) - t25 * t10 * t7) * MB1110(1) + t5)
     & + t19) * t42 - t15) + (t30)**2 * MB1111D2(0) * t42) - t24) - t31

           intHs16s25s26s34s56x1311D4eps0 = result
       end function intHs16s25s26s34s56x1311D4eps0

       function intHs16s25s26s34s56x1311D4eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1311D4eps1
           type(dd_complex) ::  t1

           type(dd_complex) :: result

      t1 = s16 + s25 + s26 - s34 + s56
      t1 = ddreal(1) / t1
      result = ddreal(1) / s16 * t1 * ((s16 + s26 - s34 + s56) * MB1110(
     &1) - ddreal(1)) / ddreal(2)

           intHs16s25s26s34s56x1311D4eps1 = result
       end function intHs16s25s26s34s56x1311D4eps1

       function intHs16s25s26s34s56x1312D6eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1312D6eps0
           type(dd_complex) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           type(dd_complex) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           type(dd_complex) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           type(dd_complex) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           type(dd_complex) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           type(dd_complex) ::  t16,t17,t18,t19,t2,t20,t21,t22,t23,t24,t25,t26
           type(dd_complex) ::  t27,t28,t29,t3,t30,t31,t32,t33,t34,t35,t36,t37
           type(dd_complex) ::  t38,t39,t4,t40,t41,t42,t43,t44,t45,t46,t47,t48
           type(dd_complex) ::  t49,t5,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59
           type(dd_complex) ::  t6,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t7
           type(dd_complex) ::  t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t8,t80
           type(dd_complex) ::  t81,t82,t83,t84,t85,t86,t87,t88,t89,t9,t90,t91
           type(dd_complex) ::  t92,t93,t94,t95,t96,t97,t98,t99

           type(dd_complex) :: result

      t1 = s16 - s25
      t2 = ddreal(9)
      t3 = ddreal(2)
      t4 = (s16)**2
      t5 = (t4)**2
      t6 = s16 * t4
      t7 = t3 * s25
      t8 = ddreal(5) * s16
      t9 = t2 * t4
      t10 = s25 * t1
      t11 = ddreal(6) * t4
      t12 = t11 - t10
      t13 = s16 + s25
      t14 = s16 + s26
      t15 = ddreal(4)
      t16 = t15 * s25
      t17 = ddreal(12) * t4
      t18 = ddreal(19) * t6
      t19 = ddreal(3) * s25
      t20 = ddreal(12) * s16
      t21 = (t20 + t19) * s25 + ddreal(5) * t4
      t22 = ddreal(20)
      t23 = ddreal(14) * t4
      t24 = t22 * t6
      t25 = t3 * s16
      t26 = ddreal(23) * t4
      t27 = -ddreal(66) * t6
      t28 = ddreal(14) * s16
      t29 = ddreal(6) * s25
      t30 = t22 * t4
      t31 = ((t29 + t28) * s25 - t30) * s25
      t32 = ddreal(36) * t6
      t33 = -t32 + t31
      t34 = ddreal(8) * s16
      t35 = ddreal(7) * t4
      t36 = (t34 + t19) * s25 - t35
      t28 = (-t7 - t28) * s25
      t37 = ddreal(40) * t6
      t38 = s25 + t25
      t39 = ddreal(27) * s16
      t40 = ddreal(7) * s25
      t41 = (s25)**2
      t42 = (t41)**2
      t43 = s25 * t42
      t44 = s25 * t41
      t45 = t15 * t41
      t46 = ddreal(24) * t4
      t47 = t46 - t45
      t48 = ddreal(3) * s16
      t49 = -s25 + t48
      t50 = ddreal(10)
      t51 = t50 * s25
      t52 = ddreal(18) * s16
      t53 = ddreal(7) * s16
      t54 = t53 + t19
      t55 = ddreal(6) * s16
      t56 = s25 - t55
      t57 = (s26)**2
      t58 = s26 * t57
      t59 = ddreal(11) * s25
      t60 = t22 * s16
      t61 = ddreal(3) * s26
      t62 = t15 * s16
      t63 = t4 * t13
      t64 = ddreal(17) * s16
      t65 = ddreal(3) * t21 * s26
      t66 = ddreal(3) * s26 * t36
      t67 = s25 * t2
      t68 = ddreal(21) * s16
      t69 = ddreal(12) * s26
      t70 = (t1)**2
      t70 = t1 * t70
      t71 = t70 * s26
      t72 = (s34)**2
      t73 = (t72)**2
      t74 = s34 * t73
      t75 = s34 * t72
      t76 = t13 * t57
      t77 = s16 * t74
      t78 = ddreal(30)
      t79 = ddreal(32) * t6
      t80 = t2 * s16
      t81 = ddreal(22) * t4
      t82 = t61 * t70
      t83 = s16 * s25
      t84 = t3 * t83 * (-t41 + t4)
      t31 = (t84 - ((-s34 * t49 - t36) * s34 - t1 * t21) * s34 + t42 - t
     &5) * s56 + s16 * ((((t20 - t16) * s25 - t81) * s25 + t18) * s25 - 
     &ddreal(6) * t5) + (((-t61 * s25 + t80 * s26 + s34 * t56 - t45 + t4
     &6) * s34 + t31 - t32 + t66) * s34 + t1 * (t41 * (t34 + t16) + ddre
     &al(24) * t63 + t65)) * s34 + t43 - t82 * t13
      t31 = t31 * s56 + t3 * (s16 * t43 - t71 * t12) + t4 * ((((s16 * t7
     &8 - t59) * s25 - ddreal(44) * t4) * s25 + t79) * s25 - t2 * t5) - 
     &ddreal(3) * t76 * t70 + ddreal(3) * t77 + (((-((t69 + t68 + t67) *
     & s16 - t7 * s26) * s34 + s16 * ((ddreal(41) * s16 + t40) * s25 + d
     &dreal(54) * t4) + s26 * (ddreal(3) * s26 * t49 + t3 * t47)) * s34 
     &+ s16 * (((t25 + t19) * s25 - t26) * s25 + t27) + s26 * (t3 * t33 
     &+ t66)) * s34 + t1 * (s16 * (((t29 + t64) * s25 - t3 * t4) * s25 +
     & ddreal(39) * t6) + (ddreal(8) * t41 * (s25 + t25) + ddreal(48) * 
     &t63 + t65) * s26)) * s34
      t32 = t70 * t58
      t34 = t14 * s34
      t63 = t34 * t1
      t65 = t71 * s16
      t66 = ddreal(16) * s16
      t85 = ddreal(15) * t6
      t86 = -t3 * (s34 * t13 + t83) + t72 + t41 + t4
      t87 = s16 + s25 - s34
      t88 = t3 * s26
      t89 = (s34 - t13 - s26) * s26 + (-t88 - t87 - s56) * s56 - t83
      t90 = t8 - t7
      t91 = -s25 + t62
      t92 = t29 + t62
      t93 = t83 * t1
      t94 = t92 * s26
      t95 = ddreal(11) * s16
      t96 = ddreal(5) * s25
      t97 = (-t96 + t95) * s25 + t46
      t98 = t8 - t19
      t99 = ddreal(3) * t4
      t100 = ddreal(24) * s16
      t101 = t3 * t58
      t102 = s25 * t14
      t103 = t15 * s26
      t104 = ddreal(13) * s16
      t105 = ddreal(3) / ddreal(2) * t92 * s26
      t106 = ddreal(13) * t4
      t107 = s25 + s26
      t108 = ddreal(8) * s26
      t109 = s16 * t42
      t110 = t71 * t91
      t111 = ddreal(17) * t4
      t112 = ddreal(3) * t70 * t57
      t113 = t4 * (((t68 - t59) * s25 - t111) * s25 + ddreal(5) * t6)
      t114 = ddreal(6) * s26
      t115 = ddreal(15) * s16
      t116 = (-t40 + t115) * s25 - t106
      t92 = t92 / ddreal(2)
      t117 = ddreal(3) * t93
      t45 = ((-t3 * s34 * (t1 * t92 + t72) + t72 * t98 - t44 + t6 - t117
     &) * s56 + t15 * (t73 + t5) + t82 + ((-(t114 + t66 + t19) * s34 + s
     &16 * (ddreal(15) * s26 + t100) + (-s26 * t2 + t95 - t96) * s25) * 
     &s34 - t1 * (s16 * (t69 + t66) + (ddreal(18) * s26 + t19 + t95) * s
     &25)) * s34 + t42 + t83 * t116) * s56 + t3 * ((-t1 * ((ddreal(11) *
     & t4 + t45) * s16 + s26 * (t105 + (t95 + t19) * s25 + ddreal(16) * 
     &t4)) + t72 * (-ddreal(16) * s16 * t14 - ddreal(3) * s26 * t107 - d
     &dreal(3) * t41 - t72 - ddreal(12) * t83)) * s34 + t109 + t110) + t
     &72 * (t72 * (t104 + t108 + t29) + (t3 * t41 + t106) * s25 + ddreal
     &(38) * t6 + s16 * t41 + (ddreal(3) * s26 * t98 + t3 * t97) * s26) 
     &+ t112 + t113
      t118 = t4 * t42
      t91 = t70 * t91 * t57
      t119 = ((t80 - t96) * s25 - t35) * s25
      t120 = s16 + s26 + s56
      t121 = t6 * t15
      t122 = ((t80 + t7) * s25 + t26) * s25
      t123 = t79 + t122
      t124 = s25 + t55
      t125 = t13 * t124
      t126 = t78 * t6
      t127 = ddreal(28) * s16
      t128 = ddreal(46) * t4
      t129 = ddreal(28) * t4
      t130 = (t39 + t29) * s25 + t129
      t68 = (t68 + t67) * s25
      t131 = s16 * t50
      t132 = (t53 + t96) * s25
      t133 = ddreal(18) * t4
      t134 = t4 * t15
      t135 = ddreal(12) * s25
      t136 = s16 * t1
      t137 = t10 * t4
      t138 = t44 - t6
      t139 = t138 * s16
      t79 = ((-t3 * s34 * (t136 * t92 + t38 * t72) - ddreal(3) * t137 - 
     &t72 * (-t125 - t72) - t139) * s56 - t3 * t74 + t4 * (s25 * t116 + 
     &t121) + ((-((ddreal(27) * s25 + t127 + t69) * s16 + ddreal(6) * s2
     &5 * t107) * s34 + t125 * t61 + t122 + t79) * s34 - t136 * (ddreal(
     &18) * s25 * s26 + t20 * s26 + t132 + ddreal(18) * t4)) * s34 + t10
     &9 + ddreal(3) * t73 * (s26 + t7 + t62) + ddreal(3) * t65) * s56 + 
     &s16 * (t113 + t112) + t3 * (s16 * (t110 + t109) + (t72 * (-t72 * (
     &t107 * t3 + t8) - (t26 - t28) * s25 - (t61 * t38 + t130) * s26 - t
     &126) - t136 * (s16 * ((t29 - t62) * s25 + t106) + s26 * (t105 + t1
     &33 + t132))) * s34) + (t72 * (t72 + (ddreal(29) * s16 + t29) * s25
     & + (t135 + t100 + t61) * s26 + ddreal(35) * t4) + t3 * t123 * s26 
     &+ s16 * ((t134 + t68) * s25 + ddreal(55) * t6) + t42 + ddreal(3) *
     & t76 * t124) * t72
      t92 = t71 * t4
      t105 = t4 * s25
      t106 = t41 + t4
      t107 = t57 + t4
      t110 = t53 + t7
      t112 = s25 + t62
      t113 = ddreal(19) * s16
      t116 = t55 + t16
      t122 = ddreal(31)
      t132 = ddreal(23) * s16
      t136 = ddreal(8) * s25
      t140 = ddreal(15) * s25
      t46 = (-t140 - t95) * s25 + t46
      t141 = -t96 + t48
      t142 = t141 * s26
      t143 = t13 * ((-t136 - t132) * s25 + t128)
      t144 = s16 * (((-t53 - t51) * s25 + t122 * t4) * s25 + ddreal(44) 
     &* t6)
      t145 = ddreal(13) * s25
      t146 = t145 - t62
      t147 = (t53 - t29) * t38
      t148 = (s16 - s26) * s26
      t149 = (ddreal(69) * s16 + ddreal(21) * s25) * s25 + ddreal(60) * 
     &t4
      t150 = t116 / ddreal(2)
      t151 = ddreal(3) * t46
      t152 = ddreal(36) * t4
      t35 = ((t141 * t72 + t3 * s34 * (-t1 * t150 + t72) - t117 - t44 + 
     &t6) * s56 + t15 * (t71 - t73 + t5) - ((-(t145 + t108 - t62) * s34 
     &- s16 * (t69 + t100) - (-s26 * t22 - t140 - t95) * s25) * s34 + t1
     & * (s16 * (ddreal(24) * s26 + t60) + (ddreal(16) * s26 + t40 + t13
     &2) * s25)) * s34 - t13 * t44 + t105 * (-t95 + t67)) * s56 - t3 * (
     &t109 - t74) + t6 * ((t140 - t113) * s25 + t35) - (((-(-t69 - t136 
     &+ t48) * s34 - (ddreal(39) * s26 + t131) * s25 + t129 + ddreal(12)
     & * t148 - ddreal(12) * t41) * s34 - ((-s16 * t122 - t136) * s25 + 
     &t26) * s25 - (ddreal(6) * t142 + t151) * s26 - ddreal(46) * t6) * 
     &s34 + t1 * (((ddreal(22) * s16 + t7) * s25 + t152) * s25 + t126 + 
     &s26 * (ddreal(6) * t116 * s26 + t149))) * s34 - t4 * t44 + ddreal(
     &6) * t70 * t57 + t82 * t112
      t1 = t35 * s56 + ddreal(3) * t70 * t112 * t57 + t3 * (t73 * (s16 *
     & (s25 + t53) + (-t136 - t114 + t48) * s26) + t65 * t110) - ddreal(
     &3) * t44 * t6 + t5 * ((t140 - t64) * s25 + t11) - (((-t72 * (-t25 
     &+ t103) - s16 * ((-t39 + t29) * s25 - t152) - ((t108 + ddreal(3) *
     & t146) * s26 - t147 * t15) * s26) * s34 - ((t142 * t15 + t151) * s
     &26 + t143 * t3) * s26 - t144) * s34 + t1 * (t15 * (s16 * t44 + (((
     &s25 + t95) * s25 + t133) * s25 + t85) * s26) + t4 * ((ddreal(25) *
     & s25 + t115) * s25 + ddreal(26) * t4) + t57 * (t108 * t150 + t149)
     &)) * s34 - t118 + t32 * t15
      t35 = s16 + s26 - s34 + s56
      t64 = ddreal(1) / (t86)**2
      t67 = ddreal(1) / (t89)**2
      t42 = t120 * (t5 * (t3 * t6 + t119) + (s26 + t25) * t72 * t73 + t7
     &9 * s56 + t42 * t6 + t72 * (-ddreal(14) * t10 * t6 + t4 * (ddreal(
     &7) * t44 + t126) + ((t125 * s26 + t123) * s26 + (((s25 + t80) * s2
     &5 + ddreal(21) * t4) * s25 + t121) * s25 + ddreal(55) * t5) * s26 
     &+ t109) + t73 * (s16 * (t4 * t78 + t68) + s26 * (t54 * (t8 + t7) +
     & t57) + ddreal(6) * t38 * t57) - t74 * (s16 * (t40 + t20) + (t88 +
     & t16 + t131) * s26) - t75 * (s16 * (((t96 + t66) * s25 + t23) * s2
     &5 + t37) + s26 * (t130 * s26 + ((t16 + t127) * s25 + t128) * s25 +
     & ddreal(60) * t6) + t101 * t38) - t63 * s16 * (s16 * ((-t80 + t40)
     & * s25 + t17) + s26 * (t94 + (s16 + t96) * s25 + t23)) + t32 * s16
     & + t92 * t90 + t91 * s16) * MB1011(1)
      t1 = s25 * (-ddreal(6) * t10 * s16 * t5 - t3 * (t138 * t5 + t74 * 
     &(t148 + t4)) + (t57)**2 * t70 + t1 * s56 + t72 * (t4 * (t41 * (-t7
     & + t131) + t24) + (((t142 + t46) * s26 + t143) * s26 + t144) * s26
     &) + t73 * (((-t103 - t136 + t48) * s26 + t25 * (s25 + t53)) * s26 
     &+ t4 * (t29 + t131)) + t75 * (t146 * t58 - t22 * t5 + t3 * t57 * (
     &-t147 + t57) - t135 * t6 - t134 * t41 - t61 * s16 * ((t80 - t7) * 
     &s25 + t17)) + t32 * t112 - t63 * (t50 * t4 * t106 - t84 + ((t116 *
     & s26 + (t40 + t113) * s25 + t23) * s26 + ((t7 + t115) * s25 + t111
     &) * s25 + ddreal(16) * t6) * s26) + t92 * t124 + s16 * t70 * t110 
     &* t57) * MB1110(1)
      t1 = -t1 - t6 * ((((-t29 + t66) * s25 - t81) * s25 + t85) * s25 - 
     &t15 * t5) - t31 * s56 - t4 * t43 - (t4 * ((t17 + t28) * s25 - t37)
     & + ((s26 * t36 + t33) * s26 + s16 * (((t25 + t19) * s25 - t26) * s
     &25 + t27)) * s26) * t72 - (t4 * ((t52 + t51) * s25 + ddreal(40) * 
     &t4) + ((s26 * t49 + t47) * s26 + s16 * t38 * (t39 + t40)) * s26) *
     & t75 - (((-t59 - t60) * s16 - t61 * t54) * s16 + t56 * t57) * t73 
     &- t77 * (t61 + t62) + t32 * t13 + t70 * t12 * t57 - t63 * (s16 * (
     &((t20 + t7) * s25 - t23) * s25 + t24) + s26 * (t21 * s26 + ((t8 + 
     &t16) * s25 + t17) * s25 + t18)) + t65 * ((-t8 + t7) * s25 + t9) + 
     &t42
      t2 = t67 * ((t19 * t87 * t72 * (s56)**3 - s25 * s34 * (t70 * (t14)
     &**2 + (((-t8 * s26 - t107 * t15 - ddreal(3) * t102 + t34) * s34 + 
     &ddreal(3) * t83 * t13 + ddreal(3) * t58 + ddreal(6) * t6 + s26 * (
     &(t80 + t40) * s26 + (t104 + t19) * s25 + t9)) * s34 - t14 * (t15 *
     & s16 * t106 - ddreal(3) * t105 + ddreal(3) * t76 + (s25 + t48) * (
     &s16 + t7) * s26 + t44)) * s34 + (-(t15 * t75 + t117 + (-(t14 * t2 
     &+ t40) * s34 + t2 * t13 * s26 + (t7 + t131) * s25 + t11) * s34 + t
     &44 - t6) * s56 - t3 * (t139 - t71) - ddreal(6) * t137 + (((-t108 -
     & t8 - t19 + s34) * s34 + (ddreal(14) * s26 + t19 + t104) * s25 + t
     &107 * t2 + t52 * s26) * s34 - t15 * (((s25 + t8) * s25 + t99) * s2
     &6 + t105) - t41 * (s25 + t55) - ddreal(7) * t6 - t76 * t2) * s34) 
     &* s56)) * MB1101(0) * t64 - t120 * (t35)**2 * MB1111D2(0))
      result = -t64 * t67 * ((s16 * I300s16(1) - MB1001(0)) * (t3 * (-t1
     &4 * t74 + t4 * t5) + t45 * s56 + t72 * (s16 * (((t53 + t7) * s25 -
     & t99) * s25 + t24) + ((s26 * t98 + t97) * s26 + t38 * ((t7 - t48) 
     &* s25 + ddreal(19) * t4)) * s26) + t73 * ((t104 + t103) * s26 + t4
     & * t50 + ddreal(6) * t102) - t75 * (s16 * ((t29 + t95) * s25 + t30
     &) + s26 * ((t66 + t19) * s26 + (t29 + t100) * s25 + ddreal(32) * t
     &4) + t101) + t32 + t118 + t91 - t63 * (t50 * t6 - ddreal(5) * t93 
     &+ s26 * (t94 + (t8 + t19) * s25 + t17)) + t65 * t90 + t119 * t6) +
     & t1 / s16) / ddreal(4) + t2 / ddreal(2)

           intHs16s25s26s34s56x1312D6eps0 = result
       end function intHs16s25s26s34s56x1312D6eps0

       function intHs16s25s26s34s56x1312D6eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1312D6eps1

           type(dd_complex) :: result

      result = ddreal(0)

           intHs16s25s26s34s56x1312D6eps1 = result
       end function intHs16s25s26s34s56x1312D6eps1

       function intHs16s25s26s34s56x1321D6eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1321D6eps0
           type(dd_complex) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           type(dd_complex) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           type(dd_complex) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           type(dd_complex) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           type(dd_complex) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           type(dd_complex) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           type(dd_complex) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           type(dd_complex) ::  t175,t176,t177,t178,t179,t18,t180,t19,t2,t20,t21,t22
           type(dd_complex) ::  t23,t24,t25,t26,t27,t28,t29,t3,t30,t31,t32,t33
           type(dd_complex) ::  t34,t35,t36,t37,t38,t39,t4,t40,t41,t42,t43,t44
           type(dd_complex) ::  t45,t46,t47,t48,t49,t5,t50,t51,t52,t53,t54,t55
           type(dd_complex) ::  t56,t57,t58,t59,t6,t60,t61,t62,t63,t64,t65,t66
           type(dd_complex) ::  t67,t68,t69,t7,t70,t71,t72,t73,t74,t75,t76,t77
           type(dd_complex) ::  t78,t79,t8,t80,t81,t82,t83,t84,t85,t86,t87,t88
           type(dd_complex) ::  t89,t9,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99

           type(dd_complex) :: result

      t1 = s16 - s25
      t2 = ddreal(2)
      t3 = t2 * s16
      t4 = ddreal(3)
      t5 = t4 * s16
      t6 = t5 - s25
      t7 = s16 + s26
      t8 = ddreal(4)
      t9 = (s16)**2
      t10 = (t9)**2
      t11 = s16 * t9
      t12 = t9 * t10
      t13 = ddreal(7) * s16
      t14 = t8 * s25
      t15 = s16 + s25
      t16 = ddreal(9) * s16
      t17 = t2 * s25
      t18 = ddreal(15) * t9
      t19 = t15 * s26
      t20 = ddreal(5)
      t21 = t20 * s16
      t22 = ddreal(10) * s16
      t23 = t22 + s25
      t24 = s16 * s25
      t25 = (s26)**2
      t26 = (t25)**2
      t27 = s26 * t26
      t28 = s26 * t25
      t29 = t2 * s26
      t30 = ddreal(15) * t11
      t31 = t8 * s16
      t32 = t20 * t9
      t33 = ddreal(11) * s16
      t34 = ddreal(17) * t9
      t35 = (-t14 + t33) * s25
      t36 = ddreal(6)
      t37 = ddreal(14) * t11
      t38 = t36 * t25
      t39 = t4 * s26
      t40 = (t1)**2
      t41 = (s25)**2
      t42 = (t41)**2
      t43 = s25 * t41
      t44 = s16 * t43
      t45 = t40 * s26
      t46 = (s34)**2
      t47 = (t46)**2
      t48 = s34 * t47
      t49 = s34 * t46
      t50 = ddreal(12) * t9
      t51 = t36 * s26
      t52 = ddreal(21)
      t53 = ddreal(15) * s16
      t54 = t4 * s25
      t55 = t9 * t52
      t56 = ddreal(9) * t11
      t57 = t45 * t36
      t58 = -t2 * (s34 * t15 + t24) + t41 + t46 + t9
      t59 = t9 * t43
      t60 = t40 * t25
      t61 = (s56)**2
      t62 = t40 * t28
      t63 = ddreal(13) * s16
      t64 = ddreal(11) * s25
      t65 = (-s34 + t15 + s26) * s26 + (t29 + s16 + s25 - s34 + s56) * s
     &56 + t24
      t66 = s16 * t15
      t67 = t2 * t41
      t68 = t36 * s16
      t69 = t13 + s25
      t70 = ddreal(16) * t11
      t71 = ddreal(14) * t9
      t72 = t20 * s25
      t73 = t5 + t72
      t74 = t73 * s26
      t75 = t11 * t4
      t76 = ddreal(34) * t11
      t77 = ddreal(12) * s25
      t78 = ddreal(36) * t9
      t79 = (ddreal(19) * s16 + t77) * s25
      t80 = t78 + t79
      t81 = ddreal(10) * t9
      t82 = ddreal(13) * t15
      t83 = ddreal(14) * s16
      t84 = ddreal(26) * t9
      t85 = (t54 + t83) * s25
      t86 = ddreal(18) * s16
      t87 = s16 * t41
      t88 = ddreal(9) * s25
      t89 = (s25 + s26) * s26 + t24
      t90 = ddreal(9) * s26
      t91 = ddreal(9) * t9
      t92 = t2 * t9
      t93 = ddreal(10) * t11
      t94 = ddreal(12) * s16
      t95 = (ddreal(33) * s25 + t94) * s25 + ddreal(51) * t9
      t96 = ddreal(8) * s25
      t97 = ddreal(15) * s25
      t98 = t8 * s26
      t99 = ddreal(22) * s25
      t100 = ddreal(36) * s16
      t101 = t10 * t1
      t102 = s16 * t42
      t103 = t9 * t41
      t104 = ddreal(8) * t28
      t105 = ddreal(7) * s25
      t106 = t36 * t41
      t107 = ddreal(20) * s26
      t108 = ddreal(17) * s16
      t109 = ddreal(12) * s26
      t110 = ddreal(7) * t9
      t111 = t9 + t41
      t73 = (((-t5 - t72 + s34) * s34 + t111 * t2 - t24 * t8) * s56 - s1
     &6 * ((-t72 + t63) * s25 - t110) - t4 * t49 + ((t98 + t82) * s34 - 
     &s16 * (t108 + t109) - (t31 + t64 + t107) * s25) * s34 + t43 + ddre
     &al(8) * t45) * s56 - s16 * ((t71 - t106) * s25 - t56) + t4 * ((t49
     & - s16 * ((-t3 + t105) * s25 + t81) - (t2 * t73 * s26 + t95 / ddre
     &al(3)) * s26 - t43) * s34 + t45 * t69) - t46 * ((t90 + t64 + t86) 
     &* s34 - ddreal(39) * t19 - t38 - t78 - t79) - t42 + ddreal(12) * t
     &60
      t73 = t73 * s56 + t101 * t20 - t2 * (-t45 * ((t31 - s25) * s25 + t
     &91) + t102) + t4 * (t47 * (t29 + t5 + s25) + t60 * t69) - (((t46 +
     & (t90 + t99 + t100) * s26 + t84 + t85) * s34 - t2 * t80 * s26 - s1
     &6 * ((t97 + t83) * s25 + ddreal(34) * t9) - t25 * (ddreal(39) * t1
     &5 + t98) - t43) * s34 + s16 * (((t21 + t96) * s25 - t92) * s25 + t
     &11 * t52) + ((t74 * t8 + t95) * s26 + t36 * (((t13 + s25) * s25 - 
     &t92) * s25 + t93)) * s26) * s34 - t103 * (t13 - t88) + t104 * t40
      t78 = t7 * t48
      t79 = t45 * s16
      t95 = s34 * t7
      t112 = t11 * t41
      t113 = t112 * t1
      t114 = t40 * t26
      t115 = s16 + s25 + s26 - s34 + s56
      t116 = t31 - s25
      t117 = t36 * t9
      t118 = s16 * (t14 + t16)
      t119 = ddreal(8) * t9
      t120 = (t21 - s25) * s25
      t121 = t36 * t11
      t122 = t60 * t4
      t123 = t36 * s25
      t124 = t9 * t8
      t125 = (t123 - t16) * s25
      t126 = t45 * t4
      t116 = t11 * ((t14 - t21) * s25 + t92) + ((-(t2 * s16 * t7 + t25) 
     &* s34 + t9 * (t17 + t68) + ((t21 + s25 + s26) * s26 + t118) * s26)
     & * s34 - t7 * (t9 * (-t54 + t68) + ((t5 - s25) * s25 + t117) * s26
     & + t87 + t2 * t15 * t25)) * s34 + ((s56 * t58 + s16 * (t124 + t125
     &) - ((-t39 - t21 - s25 + s34) * s34 + t19 * t36 + t119 + t120) * s
     &34 - t43 + t126) * s56 - t2 * ((t46 * t7 + (t19 * t4 + t119 + t120
     &) * s26 + t121) * s34 + t44 - t45 * t116) + t46 * ((t17 + t39 + t2
     &2) * s26 + t118) + t9 * ((t88 - t94) * s25 + t32) + t122) * s56 - 
     &t59 + t62 + t79 * (-t17 + t21) + t60 * t116
      t118 = ddreal(8) * s16
      t120 = ddreal(15) * t10
      t33 = t17 + t33
      t127 = -t9 + t41
      t128 = s25 * t127
      t129 = ddreal(28) * t9
      t130 = ddreal(16) * s25
      t131 = ddreal(23) * t9
      t132 = t21 + t105
      t133 = t132 * s26
      t134 = t9 * t4
      t135 = ddreal(27) * s16
      t136 = ddreal(62) * t11
      t137 = ((t17 + t135) * s25 + ddreal(32) * t9) * s25
      t138 = ddreal(19) * s25
      t139 = ddreal(63) * t9
      t140 = (ddreal(38) * s16 + t138) * s25
      t141 = t139 + t140
      t142 = ddreal(23) * s16
      t143 = ddreal(20) * s25
      t144 = t142 + t143
      t26 = t26 + t44
      t145 = ddreal(20) * t9
      t146 = ddreal(33) * s16
      t147 = ddreal(17) * s25
      t148 = ddreal(20) * s16
      t149 = (t13 - s25) * s25 + t18
      t150 = ddreal(32) * s16
      t151 = ddreal(51) * t11
      t124 = ((t72 + t150) * s25 - t124) * s25
      t152 = t3 + t14
      t153 = t152 * s25
      t154 = t110 + t153
      t155 = t11 * s25
      t156 = ddreal(24) * t9
      t157 = (t54 + t53) * s25
      t158 = t45 * t149
      t159 = t62 + t59
      t160 = ddreal(10) * s25
      t161 = t160 * t10
      t122 = t122 * t33
      t162 = ddreal(18) * s26
      t163 = (-ddreal(9) * t41 + t131) * s25
      t164 = t36 * t47
      t165 = ddreal(18) * t60
      t166 = ddreal(8) * s26
      t167 = (-t17 - t13) * s25
      t168 = ddreal(12) * t45
      t169 = t2 * s34
      t170 = t24 * t36
      t171 = t4 * t111
      t104 = ((-((t21 + t105 - t169) * s34 + t170 - t171) * s56 - (t145 
     &+ t167) * s25 - ((s34 * t36 - t142 - t143 - t166) * s34 + t8 * (dd
     &real(7) * s25 * s26 + t21 * s26 + t153 + ddreal(7) * t9)) * s34 + 
     &ddreal(11) * t11 + t168) * s56 - s16 * (-t30 + t163) - (((t138 + t
     &162 + t146) * s34 - (ddreal(69) * s16 + ddreal(60) * s25 + t109) *
     & s26 - t139 - t140) * s34 + (t133 * t36 + ddreal(12) * t154) * s26
     & + t124 + t151) * s34 - t42 + t164 + t165 + t126 * t33) * s56 - dd
     &real(9) * t11 * t127 - t2 * (t49 * (t46 + (t138 + t90 + t146) * s2
     &6 + t156 + t157) + t102 - t158) + ddreal(12) * t159 - ((-t46 * (t1
     &23 + t108 + t109) - s26 * (t4 * t144 * s26 + t2 * t141) - t136 - t
     &104 - t137) * s34 - t2 * (-(t151 + t124) * s26 + t155) + ddreal(38
     &) * t10 + ddreal(12) * t87 * t15 + ddreal(12) * t154 * t25 + t8 * 
     &t132 * t28) * s34 - t161 + t122
      t108 = s16 * t10 * s25
      t124 = t62 * t33
      t132 = t79 * ((-t17 + t118) * s25 + t91)
      t139 = t114 * t4
      t104 = t113 * t20 - t2 * (t12 - t78) - (((((ddreal(17) * s26 + t22
     &) * s16 + t36 * t89) * s34 - t25 * (t138 + t146) - t36 * ((t25 + (
     &t21 + s25) * s25 + t119) * s26 + t87) - t9 * (t147 + t148)) * s34 
     &+ t2 * t26 + t9 * ((t105 + t53) * s25 + t145) + ((s26 * t144 + t14
     &1) * s26 + t136 + t137) * s26) * s34 - t7 * (s16 * (((-t31 + t105)
     & * s25 + t134) * s25 + t93) + ((t133 + (s16 + t130) * s25 + t131) 
     &* s26 + (ddreal(16) * t41 + t129) * s16 + t128 * t20) * s26)) * s3
     &4 - t104 * s56 + t25 * ((t41 * (-t16 + s25) + ddreal(23) * t11) * 
     &s25 - t120) + t42 * t9 + t108 - t124 - t132 - t139
      t131 = s16 + s26 + s56
      t133 = (t17 + t5) * s25 + t110
      t136 = t133 * s26
      t137 = ddreal(16) * s16
      t138 = ((ddreal(39) * s16 + t88) * s25 + ddreal(48) * t9) * s25 + 
     &ddreal(102) * t11
      t140 = ddreal(26) * s16
      t141 = (t64 + t140) * s25 + ddreal(42) * t9
      t142 = t54 + t21
      t143 = ddreal(22) * t9
      t144 = ddreal(78) * t9
      t146 = (ddreal(52) * s16 + t97) * s25
      t151 = t144 + t146
      t153 = t3 + s25
      t154 = ddreal(40) * t11
      t172 = ddreal(55) * t9
      t100 = (t123 + t100) * s25
      t173 = ddreal(63) * t10
      t174 = (((t17 + t63) * s25 + t143) * s25 - t70) * s25
      t91 = (t41 * t8 + t91) * s16
      t175 = t128 + t91
      t176 = ddreal(37) * s16
      t177 = ddreal(16) * t9
      t178 = s16 * t111
      t179 = t11 + t43
      t180 = t9 * s25
      t124 = t60 * t149 - t180 * t179 + t124 + t132
      t33 = ((t178 * t4 - ((-t142 + s34) * s34 + t133) * s34 - t117 * s2
     &5) * s56 + s16 * (s16 * ((t105 - t148) * s25 + ddreal(11) * t9) + 
     &t168) + ((t4 * t46 + s16 * (ddreal(42) * s16 + t107) + (t64 + t109
     & + t140) * s25) * s34 - t8 * (t136 + t128 + t91)) * s34 + t2 * (-t
     &49 * (t29 + t72 + t22) + t44)) * s56 + s16 * (s16 * (t30 - t163) +
     & t165) - t4 * (t46 * (t49 - ((t54 + t63) * s25 + t177) * s25 - t76
     & - s26 * (t2 * t142 * s26 + t141)) - t79 * t33) - (t46 * (-(t90 + 
     &t64 + t135) * s34 + (ddreal(60) * s16 + ddreal(30) * s25 + t51) * 
     &s26 + t144 + t146) + (t136 * t36 + ddreal(12) * t175) * s26 + t173
     & + t174) * s34 - t102
      t33 = t33 * s56 + (t122 - t161) * s16 - ddreal(9) * t10 * t127 + t
     &2 * ((t158 - t102) * s16 + (t46 * (-t46 * (t17 + t39 + t13) - s16 
     &* ((t147 + t176) * s25 + ddreal(50) * t9) - (ddreal(15) * t153 * s
     &26 + t151) * s26 - t2 * (t28 + t43)) - s16 * ((t41 * (t17 + t118) 
     &- t75) * s25 + ddreal(23) * t10) - s26 * (t25 * t2 * t133 + t173 +
     & t174) - t38 * t175) * s34) + (t46 * (t46 + (ddreal(54) * s16 + t9
     &0 + t99) * s26 + t100 + t172) + t25 * t4 * t141 + s16 * (((t176 + 
     &t130) * s25 + ddreal(46) * t9) * s25 + ddreal(95) * t11) + (t36 * 
     &t138 / ddreal(3) + t142 * t8 * t25) * s26 + t42) * t46 + ddreal(12
     &) * s16 * t159
      t33 = s16 * t124 + ((((((t3 + s26) * s34 - s16 * (t105 + t94) - (t
     &39 + t14 + t83) * s26) * s34 + s16 * ((t135 + t88) * s25 + ddreal(
     &30) * t9) + t28 * t4 + s26 * ((t64 + t135) * s26 + t100 + t172)) *
     & s34 - s16 * (((t72 + t148) * s25 + ddreal(38) * t9) * s25 + t154)
     & - ((t25 + t151) * s26 + ((ddreal(34) * s16 + t14) * s25 + ddreal(
     &74) * t9) * s25 + ddreal(100) * t11) * s26 - ddreal(10) * t153 * t
     &28) * s34 + t9 * (((t105 + t118) * s25 + t143) * s25 + ddreal(30) 
     &* t11) + (((s26 * t142 + t141) * s26 + t138) * s26 + (((s25 + t137
     &) * s25 + ddreal(37) * t9) * s25 + ddreal(46) * t11) * s25 + ddrea
     &l(95) * t10) * s26 + t102) * s34 - t7 * (s16 * (((-t119 - t167) * 
     &s25 + t75) * s25 + ddreal(12) * t10) + ((t136 + ((t14 + t83) * s25
     & - t110) * s25 + ddreal(29) * t11) * s26 + ddreal(9) * t24 * t127 
     &+ t41 * (t119 + t67) + ddreal(34) * t10) * s26)) * s34 + t33 * s56
      t83 = (-t118 - t160) * s25 + t18
      t91 = t32 - t157
      t99 = t21 + t123
      t100 = ddreal(18) * s25
      t102 = -t72 - t13
      t52 = ((-s16 * t52 - t96) * s25 - t129) * s25
      t109 = t154 + t52
      t110 = (ddreal(27) * s25 + t140) * s25
      t119 = ddreal(31) * t9
      t122 = t110 + t119
      t123 = ddreal(25) * s16
      t124 = -t100 - t123
      t127 = ddreal(52) * t9
      t77 = (t137 + t77) * s25
      t128 = -t117 + t77
      t129 = t105 + t16
      t130 = t20 * t11
      t118 = ((t97 + t118) * s25 + t156) * s25 + t130
      t132 = (s16 + t88) * s25 + t50
      t133 = -ddreal(10) * t102
      t72 = (((t72 + t13 - t169) * s34 + t170 - t171) * s56 + s16 * ((t3
     &1 + t105) * s25 - t32) - t36 * (-t49 + t43) + s34 * (-(ddreal(10) 
     &* s26 + t123 + t100) * s34 + s16 * (ddreal(24) * s16 + ddreal(35) 
     &* s26) + (ddreal(25) * s26 + t100 + t3) * s25) - ddreal(15) * t45)
     & * s56 + s16 * (((-t88 + t150) * s25 - ddreal(25) * t9) * s25 + t1
     &30) - t4 * (-t49 * (t105 + t166 + t16) + t42) - t164 + s34 * (-((d
     &dreal(100) * s16 + ddreal(72) * s25 + t107) * s26 + t110 + t119) *
     & s34 + (ddreal(15) * t41 + t156) * s25 + ddreal(8) * t132 * s26 + 
     &t133 * t25 + t130 + ddreal(8) * t87) - t45 * t8 * t99 - ddreal(30)
     & * t60
      t52 = t72 * s56 + s16 * ((((t94 - t160) * s25 + t55) * s25 - ddrea
     &l(38) * t11) * s25 + t120) + t2 * t48 + (((-(t13 + t96 + t162) * s
     &34 + (ddreal(81) * s16 + ddreal(63) * s25 + ddreal(36) * s26) * s2
     &6 - t117 + t77) * s34 + ((t36 * t124 - t107) * s26 - t122 * t4) * 
     &s26 + t154 + t52) * s34 + (((ddreal(22) * s16 + t17) * s25 - t84) 
     &* s25 + ddreal(58) * t11) * s25 + ((t133 * s26 + ddreal(12) * t132
     &) * s26 + t118 * t4) * s26 - ddreal(44) * t10) * s34 - t38 * t40 *
     & t99 + t126 * t91 - ddreal(30) * t62
      t5 = (((((-t13 - t96) * s26 + t3 * (t16 + s25)) * s26 - t36 * (t28
     & - t180) + t93) * s34 + t28 * (t129 * t4 + t51) + t9 * ((-t14 - t8
     &6) * s25 - t145) - s26 * (s16 * (-t125 + t127) - s26 * t128)) * s3
     &4 + t10 * (t100 + t148) - t2 * (t103 * t15 + t27) + (s16 * (((t21 
     &- t160) * s25 - t134) * s25 + ddreal(68) * t11) + ((s26 * t124 - t
     &122) * s26 + t109) * s26) * s26) * s34 - t7 * (ddreal(10) * t9 * t
     &179 + (((t102 * s26 + (t5 - t100) * s25 - t34) * s26 + ((t22 - t97
     &) * s25 - ddreal(27) * t9) * s25 + ddreal(12) * t11) * s26 + ((t16
     &7 + t177) * s25 - ddreal(31) * t11) * s25 + ddreal(32) * t10) * s2
     &6 + t24 * ((-t67 - t177) * s25 + t121))
      t5 = t5 * s34 + (t52 * s56 - ddreal(11) * t103 * t111 + t11 * (ddr
     &eal(27) * t43 + t93) + t2 * (t47 * ((-t13 - t96) * s26 + t24 - ddr
     &eal(9) * t25 + ddreal(9) * t9) + t79 * t83) - ddreal(15) * t108 - 
     &ddreal(15) * t114 + (((t46 * (-t3 + t98) + (t106 - t127) * s16 + (
     &t2 * t128 + ddreal(24) * t25) * s26 + ddreal(9) * t129 * t25 - ddr
     &eal(9) * t180) * s34 + (t124 * t8 * t25 + t2 * t109) * s26 - t4 * 
     &(t122 * t25 + t155) + t9 * (t20 * t41 + ddreal(68) * t9) - ddreal(
     &10) * t26) * s34 + s16 * ((t43 * t8 + ddreal(25) * t11) * s25 - dd
     &real(42) * t10) + (t25 * (-t102 * t20 * s26 + ddreal(8) * t132) + 
     &(((ddreal(44) * s16 + t14) * s25 - t127) * s25 + ddreal(116) * t11
     &) * s25 - ddreal(88) * t10) * s26 - t4 * (-t118 * t25 + t59)) * s3
     &4 + t60 * (t4 * t91 - t98 * t99)) * s56 + t45 * ((s16 * t83 + (-s2
     &6 * t99 + t91) * s26) * s26 + t9 * ((t21 - t64) * s25 + t81))
      t26 = t1 * s26
      t43 = ddreal(1) / t58
      t52 = ddreal(1) / s16
      t65 = ddreal(1) / (t65)**2
      t72 = ddreal(1) / (t115)**2
      t5 = -s25 * (t2 * (t11 * t10 - t48 * ((s16 - s26) * s26 + t9)) + t
     &5 + t112 * ((-t14 + t94) * s25 - t81) - t4 * t40 * t27) * MB1110(1
     &) + t131 * (-t20 * t101 * t41 + (t12 * t2 + t139) * s16 + t33) * M
     &B1011(1)
      t5 = t104 * MB1001(0) + t5 * t52 + t8 * t113 - t73 * s56 + t25 * (
     &(t41 * (s25 - t68) + t37) * s25 - ddreal(9) * t10) - t46 * (t9 * (
     &(t17 + t16) * s25 + t81) + ((t25 + t80) * s26 + ((s25 + t53) * s25
     & + t71) * s25 + t76) * s26 + t44 + t82 * t28) - t47 * ((t21 + t90)
     & * s16 + t4 * t89) + t49 * (t4 * (t28 + t87) + t9 * (t22 + t88) + 
     &s26 * ((t64 + t86) * s26 + t84 + t85)) - t9 * (t10 - t42) + t78 - 
     &t79 * (t20 * t66 - t67) - t62 * t69 + t95 * (t20 * s16 * (t41 * (s
     &25 - s16) + t11) + (((t64 - s16) * s25 + t71 + t74) * s26 + ((t54 
     &+ t22) * s25 - t32) * s25 + t70) * s26 + t75 * s25) - t114 * t2
      t3 = t65 * t43 * ((s16 * (t9 * ((t64 - t63) * s25 + t32) + t57 * (
     &t3 - s25)) + t2 * (t58 * s56 * t61 + t62) - t4 * (-t60 * t6 + t59)
     & - ((((t29 + t31) * s26 + t32) * s34 - t8 * s16 * ((t17 + t21) * s
     &26 + t24) - (t29 + t23) * t25 - t30) * s34 + t7 * (s16 * ((t17 - t
     &16) * s25 + t18) + (t19 * t8 + (-t14 + t13) * s25 + ddreal(13) * t
     &9) * s26)) * s34 - (-(((-t54 + t53) * s25 - t55) * s25 - t2 * t49 
     &- s34 * (-(t51 + t22 + s25) * s34 + ddreal(12) * t19 + t34 + t35) 
     &+ t56 + t57) * s56 + t2 * s34 * (-((t14 + t22) * s16 + (t39 + t23)
     & * s26) * s34 + t38 * t15 + (t34 + t35) * s26 - t24 * t15 + t37) +
     & t36 * (-t45 * (s26 + t6) + t44) - t9 * ((-ddreal(30) * s16 + ddre
     &al(24) * s25) * s25 + t50) + t8 * t7 * t49) * s56) * t52 + t116 * 
     &I300s16(1) + t5 * t72)
      result = t3 / ddreal(4) + t52 * (t65 * (s25 * s34 * (-t2 * s16 * (
     &-t45 + t180) + ((s26 * t7 + t9) * s34 - t7 * ((t17 + s16) * s26 + 
     &t24 + t92)) * s34 + ((-t2 * s25 * (s16 + s34) + (s34 - s16) * s34 
     &+ t41 + t9) * s56 + t2 * (t45 + t178) - s34 * (s26 * t152 - (t29 +
     & s16) * s34 + t4 * t66) - t180 * t8) * s56 + t111 * t9 + t60) * MB
     &1101(0) * t43 + t131 * (t2 * ((-t95 + t26) * s16 + ((-s25 - s34 + 
     &s16) * s16 + t26) * s56) + (s16 * t1 + t46) * s16 + t1 * (t61 + t2
     &5)) * MB1111D2(0)) + s25 * (s16 + s26 - s34 + s56) * MB1110(0) * t
     &72) / ddreal(2)

           intHs16s25s26s34s56x1321D6eps0 = result
       end function intHs16s25s26s34s56x1321D6eps0

       function intHs16s25s26s34s56x1321D6eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1321D6eps1
           type(dd_complex) ::  t1,t2

           type(dd_complex) :: result

      t1 = s16 + s25 + s26 - s34 + s56
      t2 = ddreal(1) / ddreal(4)
      t1 = ddreal(1) / (t1)**2
      result = ddreal(1) / s16 * t1 * (t2 * (s16 - s25 + s26 - s34 + s56
     &) + s25 * (s16 + s26 - s34 + s56) * MB1110(1) / ddreal(2))

           intHs16s25s26s34s56x1321D6eps1 = result
       end function intHs16s25s26s34s56x1321D6eps1

       function intHs16s25s26s34s56x1411D6eps0()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1411D6eps0
           type(dd_complex) ::  t1,t10,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109
           type(dd_complex) ::  t11,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t12
           type(dd_complex) ::  t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130
           type(dd_complex) ::  t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141
           type(dd_complex) ::  t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152
           type(dd_complex) ::  t153,t154,t155,t156,t157,t158,t159,t16,t160,t161,t162,t163
           type(dd_complex) ::  t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174
           type(dd_complex) ::  t175,t176,t177,t178,t179,t18,t180,t181,t182,t183,t184,t185
           type(dd_complex) ::  t186,t187,t188,t189,t19,t190,t191,t192,t193,t194,t195,t196
           type(dd_complex) ::  t197,t198,t199,t2,t20,t200,t201,t202,t203,t204,t205,t206
           type(dd_complex) ::  t207,t208,t209,t21,t210,t211,t212,t213,t214,t215,t216,t217
           type(dd_complex) ::  t218,t219,t22,t220,t221,t222,t223,t224,t225,t226,t227,t228
           type(dd_complex) ::  t229,t23,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239
           type(dd_complex) ::  t24,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249,t25
           type(dd_complex) ::  t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t26,t260
           type(dd_complex) ::  t261,t262,t263,t264,t265,t27,t28,t29,t3,t30,t31,t32
           type(dd_complex) ::  t33,t34,t35,t36,t37,t38,t39,t4,t40,t41,t42,t43
           type(dd_complex) ::  t44,t45,t46,t47,t48,t49,t5,t50,t51,t52,t53,t54
           type(dd_complex) ::  t55,t56,t57,t58,t59,t6,t60,t61,t62,t63,t64,t65
           type(dd_complex) ::  t66,t67,t68,t69,t7,t70,t71,t72,t73,t74,t75,t76
           type(dd_complex) ::  t77,t78,t79,t8,t80,t81,t82,t83,t84,t85,t86,t87
           type(dd_complex) ::  t88,t89,t9,t90,t91,t92,t93,t94,t95,t96,t97,t98
           type(dd_complex) ::  t99

           type(dd_complex) :: result

      t1 = s16 - s25
      t2 = ddreal(15) * s16
      t3 = ddreal(8) * s25
      t4 = t3 + t2
      t5 = ddreal(4)
      t6 = (s16)**2
      t7 = (t6)**2
      t8 = (t7)**2
      t9 = s16 * t6
      t10 = t9 * t7
      t11 = t6 * t7
      t12 = s16 * t7
      t13 = t5 * s25
      t14 = ddreal(29) * s16
      t15 = ddreal(43) * t6
      t16 = s16 + s25
      t17 = ddreal(2)
      t18 = t17 * s25
      t19 = ddreal(46) * t6
      t20 = t5 * s16
      t21 = s25 + t20
      t22 = ddreal(11) * s25
      t23 = ddreal(13) * s16
      t24 = t22 + t23
      t25 = s16 + s26
      t26 = ddreal(9)
      t27 = ddreal(43) * s16
      t28 = t26 * s25
      t29 = ddreal(30) * t6
      t30 = ddreal(73) * s25
      t31 = ddreal(141) * t9
      t32 = ddreal(98) * s16
      t33 = ddreal(82) * s25
      t34 = ddreal(96) * t6
      t35 = ddreal(23) * s16
      t36 = ddreal(27) * s25
      t37 = t36 + t35
      t38 = (s25)**2
      t39 = (t38)**2
      t40 = s25 * t38
      t41 = t40 * t39
      t42 = t38 * t39
      t43 = s25 * t39
      t44 = t6 + t38
      t45 = t6 * s25
      t46 = t45 * t44
      t47 = (s26)**2
      t48 = (t47)**2
      t49 = s26 * t48
      t50 = s26 * t47
      t51 = ddreal(52) * t6
      t52 = t17 * s26
      t53 = ddreal(119) * s16
      t54 = ddreal(6)
      t55 = t54 * s16
      t56 = ddreal(39) * s25
      t57 = ddreal(19) * t6
      t58 = (t56 + t55) * s25 + t57
      t59 = t54 * t6
      t60 = ((ddreal(69) * s16 + ddreal(168) * s25) * s25 + t59) * s25
      t61 = ddreal(157) * t9
      t62 = ddreal(5) * s26
      t63 = t62 * t58
      t64 = ddreal(172) * s25
      t65 = ddreal(45) * s26
      t66 = ddreal(30) * s26
      t67 = ddreal(45) * s16
      t68 = ddreal(88) * s25
      t69 = (t1)**2
      t70 = (t69)**2
      t71 = t1 * t69
      t72 = ddreal(34)
      t73 = s16 * t72
      t74 = ddreal(16) * s25
      t75 = (s34)**2
      t76 = (t75)**2
      t77 = s34 * t75
      t78 = t77 * t76
      t79 = t75 * t76
      t80 = s34 * t76
      t81 = ddreal(24) * t6
      t82 = (-t73 + t74) * s25
      t83 = ddreal(35) * s26
      t84 = ddreal(25) * s25
      t85 = t26 * s16
      t86 = -t84 - t85
      t87 = s16 * s25
      t88 = t87 * t44
      t89 = t7 + t39
      t90 = ddreal(42) * t6
      t91 = ddreal(762) * s16
      t92 = ((ddreal(219) * s25 + t91) * s25 + ddreal(846) * t6) * s25 +
     & ddreal(711) * t9
      t93 = (ddreal(750) * s16 + ddreal(492) * s25) * s25 + ddreal(714) 
     &* t6
      t94 = ddreal(10) * t37 * s26
      t95 = ddreal(18) * s25
      t96 = ddreal(49)
      t97 = t96 * t9
      t98 = (((ddreal(309) * s16 + ddreal(200) * s25) * s25 + t90) * s25
     & + t97) * s25
      t99 = ddreal(408) * t7
      t100 = t98 + t99
      t101 = t60 + t61
      t102 = t9 * t38
      t103 = ddreal(3)
      t104 = ddreal(66) * s25
      t105 = (ddreal(361) * s16 + t104) * s25 + ddreal(191) * t6
      t106 = t100 * t103
      t107 = ddreal(10) * t58
      t108 = ddreal(177) * t7
      t109 = ddreal(475) * t12
      t110 = ddreal(340) * t9
      t111 = ((ddreal(409) * s16 + ddreal(290) * s25) * s25 + ddreal(247
     &) * t6) * s25
      t112 = t111 + t110
      t113 = ddreal(95) * t6
      t114 = (t64 + t53) * s25 + t113
      t115 = -ddreal(10) * t86 * s26
      t116 = t112 * t103
      t117 = (ddreal(774) * s16 + ddreal(705) * s25) * s25 + ddreal(486)
     & * t6
      t118 = t68 + t67
      t119 = ddreal(60)
      t120 = s26 * t119
      t121 = ddreal(303) * s25
      t122 = ddreal(189) * s16
      t123 = ddreal(108) * s26
      t124 = ddreal(102) * s25
      t125 = ddreal(54)
      t126 = ddreal(27) * s16
      t127 = ddreal(38) * s25
      t128 = t125 * s26
      t129 = ddreal(12)
      t130 = ddreal(37) * t6
      t131 = t129 * t70
      t132 = ddreal(70) * t70
      t133 = t103 * s26
      t134 = t54 * t78
      t135 = ddreal(2) / ddreal(3) * t93
      t136 = ddreal(254) * s16
      t137 = t101 * t5
      t138 = t5 * t114
      t139 = ddreal(90)
      t140 = ddreal(110) * t6
      t141 = t70 * s26
      t142 = t129 * s25
      t143 = ddreal(11) * s16
      t144 = t5 * t9
      t145 = ddreal(36)
      t146 = t145 * s25
      t147 = t9 * t54
      t148 = ((s25 * t105 - t147) * s25 + t108) * s25 + t109
      t149 = ((((ddreal(201) * s16 + ddreal(100) * s25) * s25 + ddreal(3
     &8) * t6) * s25 + ddreal(18) * t9) * s25 + ddreal(162) * t7) * s25 
     &+ ddreal(265) * t12
      t150 = ddreal(260) * s25
      t151 = (((ddreal(988) * s16 + t150) * s25 + ddreal(770) * t6) * s2
     &5 + ddreal(778) * t9) * s25 + ddreal(992) * t7
      t91 = ((ddreal(300) * s25 + t91) * s25 + ddreal(620) * t6) * s25 +
     & ddreal(568) * t9
      t152 = ddreal(215) * t9
      t153 = ((ddreal(199) * s16 + ddreal(140) * s25) * s25 + ddreal(235
     &) * t6) * s25 + t152
      t154 = ddreal(204) * s25
      t155 = (ddreal(314) * s16 + t154) * s25 + ddreal(194) * t6
      t156 = t121 + t122
      t157 = ddreal(64) * s25
      t158 = (ddreal(75) * s16 + t157) * s25 + ddreal(76) * t6
      t159 = s16 * t129
      t160 = ddreal(11) * t6
      t161 = (((((-(t125 * s26 * t25 + (ddreal(76) * s26 + t159) * s25 +
     & t160) * s34 + s16 * t158 + ddreal(72) * t50 + s26 * (s26 * t156 +
     & t155)) * s34 - s16 * t153 - t50 * (t118 * t5 + t66) - s26 * (s26 
     &* t117 + t91)) * s34 + s16 * ((((ddreal(270) * s16 + ddreal(160) *
     & s25) * s25 + ddreal(204) * t6) * s25 + ddreal(334) * t9) * s25 + 
     &ddreal(320) * t7) + t47 * (-ddreal(5) * t86 * t47 + t116) + s26 * 
     &(t138 * t47 + t151)) * s34 - s16 * t149 - (((t137 + t63) * s26 + t
     &106) * s26 + t148 * t17) * s26) * s34 + t69 * (s16 * ((((ddreal(14
     &3) * s16 + ddreal(32) * s25) * s25 + ddreal(198) * t6) * s25 + ddr
     &eal(191) * t9) * s25 + ddreal(116) * t7) + ((s26 * t135 + t92) * s
     &26 + (((ddreal(354) * s16 + t146) * s25 + ddreal(636) * t6) * s25 
     &+ ddreal(672) * t9) * s25 + ddreal(466) * t7) * s26 + ddreal(5) * 
     &t37 * t48)) * s34
      t162 = t70 * t16 * ((t14 + t18) * s25 + t19)
      t163 = t103 * t70
      t164 = t10 * s25
      t165 = t9 * s25
      t166 = -ddreal(86)
      t167 = ddreal(68)
      t168 = ddreal(44)
      t169 = ddreal(80) * t7
      t170 = ddreal(41) * s16
      t171 = ddreal(26) * s25
      t172 = ddreal(69) * t6
      t173 = t54 * t49
      t91 = t6 * (((t170 + t171) * s25 + t172) * s25 + t119 * t9) + (s16
     & * t153 + ((s26 * t118 + t117 / ddreal(3)) * s26 + t91 / ddreal(2)
     &) * s26) * s26 + t173
      t153 = t54 * s25
      t174 = t70 * t49
      t175 = t8 * s25
      t60 = (t54 * t76 - ddreal(28) * t88 + ddreal(7) * t89 - ((-s34 * t
     &86 - t58) * s34 + t37 * t69) * s34 + t90 * t38) * s56 + (((t82 - t
     &81) * s25 + ddreal(116) * t9) * s25 - ddreal(104) * t7) * s25 - dd
     &real(18) * t80 + ddreal(30) * t12 - (((-(t68 + t67 + t66) * s34 + 
     &s16 * (ddreal(95) * s16 + t65) + (ddreal(125) * s26 + t53 + t64) *
     & s25) * s34 - t60 - t61 - t63) * s34 + t69 * (s16 * (ddreal(115) *
     & s26 + t53) + (ddreal(125) * s16 + ddreal(135) * s26 + t33) * s25)
     &) * s34 + t83 * t70
      t23 = t60 * s56 + (((((t22 + t23) * s25 - t140) * s25 + t139 * t9)
     & * s25 + ddreal(95) * t7) * s25 - ddreal(151) * t12) * s25 + ddrea
     &l(52) * t11 + ddreal(18) * t79 - (((((ddreal(63) * s16 + ddreal(10
     &1) * s25 + ddreal(72) * s26) * s34 - (ddreal(258) * s16 + ddreal(2
     &35) * s25) * s25 - (ddreal(180) * s16 + ddreal(352) * s25 + t120) 
     &* s26 - ddreal(162) * t6) * s34 + s26 * (t138 + t115) + t111 + t11
     &0) * s34 - (t107 * s26 + t137) * s26 - t98 - t99) * s34 + t69 * ((
     &(t136 + t30) * s25 + ddreal(282) * t6) * s25 + ddreal(237) * t9 + 
     &s26 * (t135 + t94))) * s34 + ddreal(8) * t141 * t4 + t132 * t47
      t23 = t23 * s56 + ((((((t35 + t18) * s25 - t130) * s25 - ddreal(76
     &) * t9) * s25 + ddreal(144) * t7) * s25 + ddreal(7) * t12) * s25 -
     & ddreal(109) * t11) * s25 + ddreal(46) * t10 - (((((-(t126 + t127 
     &+ t128) * s34 + (ddreal(157) * s16 + t124) * s25 + (t123 + t121 + 
     &t122) * s26 + ddreal(97) * t6) * s34 - ((ddreal(381) * s16 + ddrea
     &l(150) * s25) * s25 + ddreal(310) * t6) * s25 - ((t118 * t54 + t12
     &0) * s26 + t117) * s26 - ddreal(284) * t9) * s34 + (((ddreal(494) 
     &* s16 + ddreal(130) * s25) * s25 + ddreal(385) * t6) * s25 + ddrea
     &l(389) * t9) * s25 + ((t114 * t54 + t115) * s26 + t116) * s26 + dd
     &real(496) * t7) * s34 - (t105 * t38 + t108) * s25 - (t107 * t47 + 
     &t106) * s26 + t54 * (-t101 * t47 + t102) - t109) * s34 + t69 * (((
     &(ddreal(177) * s16 + t95) * s25 + ddreal(318) * t6) * s25 + ddreal
     &(336) * t9) * s25 + ddreal(233) * t7 + ((t93 + t94) * s26 + t92) *
     & s26)) * s34 + t131 * t4 * t47 + t132 * t50 + t133 * t70 * t21 * t
     &24 - t134
      t23 = t23 * s56 + (-t129 * t78 + t162 * t17) * s26 - t161 + ddreal
     &(21) * t8 + ddreal(8) * t70 * t4 * t50 + ddreal(35) * t70 * t48 + 
     &t163 * t21 * t24 * t47 + t20 * t41 - ddreal(41) * t164 + ddreal(62
     &) * t12 * t40 - t97 * t43 - ddreal(17) * t11 * t38 + ddreal(7) * t
     &7 * t39 + ddreal(13) * t6 * t42
      t60 = t9 * t42
      t61 = t69 * t25
      t30 = t61 * (s16 * (t38 * (ddreal(14) * t38 + t51) + ddreal(24) * 
     &t7) + ddreal(39) * t46 + ((t37 * s26 + (t33 + t32) * s25 + t34) * 
     &s26 + ((ddreal(172) * s16 + t30) * s25 + ddreal(184) * t6) * s25 +
     & t31) * t47 + t52 * t16 * (((t28 + t27) * s25 + t29) * s25 + ddrea
     &l(46) * t9)) * s34
      t33 = t141 * s16
      t37 = t162 * t47
      t14 = t23 * s56 + t75 * (t6 * (t119 * t12 + t40 * (t145 * t38 - t5
     &9)) + t72 * t165 * ((t6 + t38) * s25 + t9) + (s16 * t149 + (((s26 
     &* t58 + t101) * s26 + t100) * s26 + t148) * s26) * s26) + t76 * t9
     &1 + t77 * (s16 * (s16 * ((((-s25 * t168 - ddreal(50) * s16) * s25 
     &- t167 * t6) * s25 + t166 * t9) * s25 - t169) - t52 * ((((ddreal(1
     &35) * s16 + ddreal(80) * s25) * s25 + ddreal(102) * t6) * s25 + dd
     &real(167) * t9) * s25 + ddreal(160) * t7)) - (((-t86 * s26 + t114)
     & * s26 + t112) * s26 + t151 / ddreal(2)) * t47) + t79 * (ddreal(18
     &) * t50 + s26 * (s16 * (t143 + t142) + s26 * (t126 + t127)) + t144
     &) - t80 * (t6 * ((ddreal(19) * s16 + t153) * s25 + t81) + ddreal(1
     &8) * t48 + (s16 * t158 + (s26 * t156 / ddreal(3) + t155 / ddreal(2
     &)) * s26) * s26) + t60 - t30 + t33 * (((t14 + t13) * s25 + t15) * 
     &s25 + ddreal(21) * t9) + t70 * t21 * t24 * t50 + t37
      t21 = ddreal(15) * s25
      t23 = ddreal(14) * s16
      t24 = t103 * t9
      t30 = t7 * t38
      t37 = s34 * t16
      t58 = -t17 * (t37 + t87) + t6 + t38 + t75
      t63 = s16 + s25 + s26 - s34 + s56
      t64 = (-s34 + t16 + s26) * s26 + (t52 + s16 + s25 - s34 + s56) * s
     &56 + t87
      t67 = ddreal(7) * s16
      t68 = ddreal(8) * s16
      t86 = s25 + t68
      t91 = ddreal(7) * s25
      t92 = ddreal(16) * s16
      t93 = t85 + t142
      t94 = ddreal(5) * s16
      t98 = ddreal(24) * s25
      t99 = t17 * t6
      t100 = ddreal(22)
      t101 = s16 * t100
      t105 = t38 * (t101 + t95) + ddreal(32) * t9
      t106 = t103 * t38
      t107 = t6 + t106
      t109 = s16 * (((t98 + t67) * s25 + t99) * s25 + ddreal(39) * t9)
      t110 = ddreal(7) * t38
      t111 = t129 * t9
      t112 = ddreal(20) * s25
      t114 = ddreal(18) * s16
      t115 = t100 * s25
      t116 = (t115 + t114) * s25 + ddreal(15) * t6
      t117 = s16 + t13
      t118 = t6 * t38
      t121 = s16 * ((t112 + t159) * s25 + ddreal(25) * t6)
      t122 = ddreal(13) * s25
      t132 = s16 * t16
      t135 = t17 * s16
      t137 = s25 + t135
      t138 = t93 / ddreal(3)
      t148 = ddreal(24) * s16
      t149 = ddreal(27) * t6
      t151 = t26 * s26
      t155 = t6 * t39
      t156 = ddreal(26) * t6
      t158 = t141 * t26
      t161 = t80 * t54
      t162 = ddreal(27) * s26
      t176 = t145 * s26
      t177 = t6 * t168
      t178 = t59 * t38
      t179 = (s56)**2
      t180 = s56 * (t179)**2
      t181 = t40 * (s16 + s25) - t7
      t82 = (s16 * ((((-ddreal(26) * s16 + t13) * s25 + t177) * s25 - dd
     &real(31) * t9) * s25 + ddreal(8) * t7) + t75 * (-t103 * t77 + t17 
     &* (t26 * (s26 * t6 + t40) + t38 * (t162 + t143) + ddreal(16) * t9)
     &) + t158 + (-t69 * ((t146 + t126) * s26 + t110 + ddreal(28) * t132
     &) + t75 * ((t122 + t55 + t151) * s34 - s16 * (t151 + t2) - (t115 +
     & t114 + t176) * s25)) * s34 + t43) * s56 + (t158 - t161) * s26 + t
     &17 * (s16 * t43 + t141 * t86) + t9 * ((-t82 - t156) * s25 + ddreal
     &(7) * t9) - (((-((t159 + t171 + t151) * s26 + t132 * t54) * s34 + 
     &(t17 * t116 + t151 * t117) * s26 + t121) * s34 - (ddreal(18) * t10
     &7 * s26 + t17 * t105) * s26 - t109) * s34 + t69 * (s16 * ((t148 + 
     &t142) * s25 + t149) + (ddreal(14) * (t137)**2 + t151 * t138) * s26
     &)) * s34 - t155
      t158 = t6 * t43
      t171 = t12 * s25
      t182 = (t92 + t18) * s25 + t160
      t183 = ddreal(17) * s16
      t184 = (t183 + s25) * s25 + ddreal(21) * t6
      t185 = t183 + t153
      t186 = ddreal(31) * s25
      t187 = ddreal(10) * s25
      t188 = ddreal(10) * t6
      t189 = ddreal(78) * t6
      t190 = -ddreal(105) * t9
      t191 = (((s16 + ddreal(29) * s25) * s25 - t189) * s25 + t190) * s2
     &5
      t192 = t108 + t191
      t193 = ddreal(51)
      t194 = ddreal(81) * s16
      t195 = ddreal(30) * s25
      t196 = t193 * t9
      t197 = ((t195 - t194) * s25 - ddreal(88) * t6) * s25
      t198 = t196 + t197
      t199 = ddreal(21) * s16
      t200 = ddreal(14) * t6
      t201 = (t28 + t199) * s25
      t202 = t200 + t201
      t203 = s25 * t1
      t204 = ddreal(16) * t6
      t205 = ddreal(93)
      t206 = ddreal(138) * t9
      t207 = ((s16 * t205 - t74) * s25 + ddreal(92) * t6) * s25
      t208 = -t206 + t207
      t209 = (ddreal(134) * s16 + t146) * s25
      t210 = t209 + t130
      t211 = ddreal(57) * s16
      t212 = ddreal(52) * s25
      t213 = t211 + t212
      t214 = ddreal(48)
      t215 = t52 * s16
      t216 = t214 * s25
      t217 = ddreal(23) * t6
      t218 = ddreal(21) * s25
      t166 = (s16 * t166 - t218) * s25
      t219 = t149 + t166
      t220 = ddreal(72) * s16
      t221 = -s25 * t167 - t220
      t222 = s16 * t38
      t223 = ddreal(30) * s16
      t224 = ddreal(41) * t6
      t225 = ddreal(33) * t9
      t226 = ((t101 + t13) * s25 + ddreal(20) * t6) * s25
      t227 = t119 * t6
      t228 = (s16 * t145 + t186) * s25
      t229 = t227 + t228
      t230 = ddreal(28) * s16
      t231 = t216 + t230
      t232 = ddreal(15) * s26
      t233 = ddreal(32) * s16
      t234 = ddreal(40) * t6
      t235 = ddreal(24) * t9
      t236 = ddreal(8) * t202
      t237 = t16 * s26
      t238 = t198 * t103
      t239 = ddreal(32) * t6
      t240 = ddreal(170) * t9
      t241 = t210 * t103
      t242 = ddreal(42) * s16
      t243 = ddreal(81) * s26
      t244 = s16 * t42
      t245 = t141 * t184 + t244
      t144 = t6 * (((((-t170 + t3) * s25 + t177) * s25 + t144) * s25 - d
     &dreal(28) * t7) * s25 + ddreal(11) * t12)
      t177 = ddreal(20) * t70 * t50
      t246 = t163 * t185 * t47
      t247 = ddreal(120) * t237
      t248 = t139 * s26
      t249 = -ddreal(216)
      t250 = ddreal(28) * s25
      t251 = t185 * s26
      t252 = s16 * (((((-t170 + t122) * s25 + t200) * s25 + ddreal(59) *
     & t9) * s25 - ddreal(67) * t7) * s25 + ddreal(21) * t12)
      t253 = t26 * t79
      t254 = ddreal(30) * t70
      t255 = t254 * t47
      t256 = ddreal(17) * s25
      t257 = ddreal(78) * t9
      t258 = ddreal(20) * t141
      t259 = ddreal(17) * t12
      t260 = ((((-t67 + t153) * s25 - t239) * s25 + t257) * s25 - ddreal
     &(62) * t7) * s25
      t261 = t87 * t17
      t29 = ((t75 * (-t17 * (ddreal(30) * t237 + t200 + t201) - ddreal(2
     &7) * t75) - ddreal(20) * t88 + ddreal(5) * t89 - s34 * (t69 * (t67
     & + t142 - t232) - t75 * (t211 + t212 + t65)) + t29 * t38) * s56 + 
     &t76 * (-t5 * (t256 + t114 + t162) + ddreal(27) * s34) - ((-((ddrea
     &l(228) * s16 + ddreal(208) * s25 + t248) * s26 + t130 + t209) * s3
     &4 - (-t236 - t247) * s26 - t196 - t197) * s34 + t69 * ((t216 - t66
     & + t230) * s26 + t227 + t228)) * s34 + t258 + t259 + t260) * s56 +
     & t103 * t69 * (-(-ddreal(10) * t50 + s26 * (t231 * s26 / ddreal(2)
     & + t229) + t225 + t226) * s34 + t251 * t69) + t75 * ((((t199 + t24
     &3 + t250) * s34 + (s16 * t249 - ddreal(162) * s26 - t154) * s26 + 
     &t149 + t166) * s34 + ((t213 * t54 + t248) * s26 + t241) * s26 - t2
     &06 + t207) * s34 + ((-t129 * t202 - t247) * s26 + t238) * s26 + t1
     &08 + t191) + t42 + t252 - t253 + t255
      t27 = t29 * s56 + t17 * t245 - ((((-((ddreal(56) * s25 + t242 + t2
     &43) * s26 - t224) * s34 - ((t129 * t221 / ddreal(4) - t123) * s26 
     &+ t17 * t219) * s26 - ddreal(20) * t87 * t1 - ddreal(115) * t9) * 
     &s34 - s16 * (((s16 * t168 - t18) * s25 - t239) * s25 - t240) - (((
     &t213 * t5 + t65) * s26 + t241) * s26 + t17 * t208) * s26) * s34 - 
     &s16 * ((((t146 - t233) * s25 - t234) * s25 - t235) * s25 + ddreal(
     &140) * t7) - (((-t237 * t119 - t236) * s26 + t238) * s26 + t17 * t
     &192) * s26) * s34 + t69 * (s16 * (((t115 + t27) * s25 + ddreal(58)
     & * t6) * s25 + ddreal(61) * t9) + (((t231 - t232) * s26 + t103 * t
     &229) * s26 + t54 * (t225 + t226)) * s26)) * s34 + t144 + t177 + t2
     &46 + t54 * (s16 - t133) * t79
      t29 = t61 * (s16 * (((t68 + t187) * s25 + t100 * t6) * s25 + t111)
     & + t50 * (t17 * (t94 + t153) - t133) + s26 * (((t186 + t148) * s25
     & + ddreal(50) * t6) * s26 + ((ddreal(35) * s16 + t142) * s25 + t14
     &5 * t6) * s25 + t97)) * s34
      t42 = ddreal(7) * t6
      t65 = (t114 - t142) * s25 - t42
      t97 = t17 * t7
      t108 = t38 * t65 + t97
      t114 = ddreal(5) * t70 * t48
      t27 = t103 * (s34 * (t103 * t75 - t37 * t5 - t261 + t44) * t180 - 
     &t164) + t7 * t108 + t27 * s56 + t75 * (t100 * t203 * t7 + ddreal(3
     &0) * t11 + t47 * (s26 * t198 + t192) + t20 * ((((t28 - t68) * s25 
     &- t188) * s25 - t147) * s25 + ddreal(35) * t7) * s26 - t129 * t16 
     &* t49 - t17 * t202 * t48 + ddreal(8) * t155 - t147 * t40) + t76 * 
     &(t9 * (ddreal(37) * s25 + t223) - ddreal(27) * t48 + ddreal(5) * s
     &16 * ((t203 * t5 + t217) * s26 + t222) + t47 * (s26 * t221 + t219)
     &) + t77 * (-ddreal(40) * t12 + t26 * t49 + ((s26 * t213 + t210) * 
     &s26 + t208) * t47 - t215 * (((-t101 + s25) * s25 + t204) * s25 + d
     &dreal(85) * t9) - t216 * t7 + t5 * t6 * t40) + t79 * ((t55 - t151)
     & * s26 + t99) + t80 * (((ddreal(7) * t138 + t162) * s26 - t224) * 
     &s26 - t6 * (t159 + t187)) - t29 + t33 * t182 + t70 * t185 * t50 + 
     &t70 * t184 * t47 + t158 * t16 + t114
      t28 = s16 + s26 + s56
      t29 = ddreal(57) * t9
      t37 = ddreal(18) * t6
      t44 = s16 + t18
      t99 = t44 * s26
      t100 = -s16 - t153
      t101 = t6 * t214
      t108 = ddreal(42) * t9
      t115 = t103 * s25
      t2 = ((((-t115 - t2) * s25 - t51) * s25 - t257) * s25 - ddreal(105
     &) * t7) * s25 + ddreal(277) * t12
      t51 = ddreal(94)
      t123 = t6 * t139
      t130 = t51 * t9
      t146 = (((-t256 - t230) * s25 - t123) * s25 - t130) * s25
      t147 = ddreal(141) * t7
      t151 = t147 + t146
      t154 = ddreal(45) * t6
      t164 = ((-t218 - t73) * s25 - t154) * s25
      t166 = t111 + t164
      t191 = t94 + t91
      t192 = t237 * t191
      t196 = s16 * ((((t100 * s25 - t101) * s25 - t108) * s25 - ddreal(1
     &8) * t7) * s25 + ddreal(195) * t12)
      t197 = ddreal(42) * t12
      t198 = ddreal(19) * s25
      t73 = (-t73 - t198) * s25
      t200 = ddreal(148) * t9
      t201 = ddreal(295) * t7
      t202 = ddreal(109) * t6
      t206 = -ddreal(335) * t7
      t168 = (((t220 + t21) * s25 + t202) * s25 - t168 * t9) * s25
      t207 = t206 + t168
      t53 = ((t53 + t216) * s25 + t123) * s25
      t123 = ddreal(104) * t9
      t208 = -t123 + t53
      t167 = s16 * t167
      t209 = (t56 + t167) * s25 + t204
      t210 = t67 + t3
      t212 = ddreal(80) * t9
      t149 = ((-t195 - t194) * s25 + t149) * s25 + t152
      t152 = (-ddreal(100) * s16 - ddreal(62) * s25) * s25 + t160
      t194 = ddreal(25) * s16
      t213 = -t186 - t194
      t74 = s16 * (((t223 - t74) * s25 + ddreal(192) * t6) * s25 + ddrea
     &l(265) * t9)
      t219 = ddreal(30) * s25 * t16
      t221 = ddreal(61) * t6
      t224 = -t221 + t219
      t226 = t148 + t127
      t227 = ddreal(42) * t7
      t228 = t37 * t38
      t170 = s16 * (ddreal(14) * s25 + t170)
      t22 = -t26 * t50 + t6 * (t22 + t23) + s26 * (-s26 * (s16 + t21) + 
     &t170)
      t34 = ((s16 * t193 - t153) * s25 + t34) * s25 + ddreal(240) * t9
      t229 = (-t13 + t55) * s25 + ddreal(17) * t6
      t230 = t229 * t5
      t231 = t166 * t5
      t191 = t16 * t191
      t236 = t191 * t62
      t237 = t103 * t208
      t238 = t209 * t5
      t241 = t103 * t152
      t243 = t54 * s26
      t247 = ddreal(64) * s16
      t248 = t44 * t47
      t249 = ddreal(52) * s16
      t257 = ddreal(10) * t192
      t262 = ddreal(10) * t210
      t263 = -ddreal(124)
      t264 = ddreal(40) * s26
      t265 = ddreal(30) * t9
      t46 = (t26 * t80 - ddreal(20) * t46 + (((-(t186 + t194 + t232) * s
     &34 + s16 * (t92 + t83) + (t56 + t167 + t264) * s25) * s34 + t111 +
     & t164 - t236) * s34 - t69 * (s16 * (t183 - t62) + (-ddreal(10) * s
     &26 - t13 + t55) * s25)) * s34 + ddreal(5) * s16 * t89 + t265 * t38
     &) * s56 + s16 * (t258 + t259 + t260) - t253 + (((((t148 + t176 + t
     &127) * s34 + (s25 * t263 - t66) * s26 + t160 - ddreal(62) * t38 - 
     &ddreal(100) * (s25 + s26) * s16) * s34 + (s26 * t262 + t238) * s26
     & - t123 + t53) * s34 + (-t257 + t231) * s26 + t146 + t147) * s34 -
     & t69 * (((t183 - t18) * s25 + t239) * s25 + (t230 - ddreal(10) * t
     &99) * s26 + t212)) * s34
      t46 = t46 * s56 + s16 * (t252 + t255) + t103 * (t33 * t185 + t78) 
     &+ (((((-(s16 + t162 + t21) * s34 + (ddreal(114) * s25 + t128 + t22
     &0) * s26 + t219 - t221) * s34 + s16 * ((-ddreal(81) * s25 + t126) 
     &* s25 + ddreal(215) * t6) + (t243 * t213 + t241) * s26 - ddreal(30
     &) * t40 - ddreal(30) * t50) * s34 + t206 + t168 + t54 * t209 * t47
     & + s26 * (t262 * t47 + t237)) * s34 + s16 * ((((-t249 - t21) * s25
     & - t189) * s25 + t190) * s25 + ddreal(277) * t7) + t103 * (s26 * t
     &151 - t43) + t47 * (t166 * t54 - t257)) * s34 - t69 * ((-ddreal(10
     &) * t248 + t34) * s26 + t54 * (s16 * t40 + t229 * t47) + t6 * ((t2
     &47 + t216) * s25 + ddreal(119) * t6))) * s34 + t244
      t34 = t46 * s56 + s16 * (t144 + t177 + t246) + t17 * s16 * t245 - 
     &((((((-(-t94 + t243) * s34 - (-t195 - t162 - t135) * s26 - t170) *
     & s34 - s16 * ((-s16 * t139 - t153) * s25 - ddreal(141) * t6) - t14
     &5 * t50 - s26 * (t103 * t226 * s26 + t17 * t224)) * s34 - t50 * (t
     &213 * t5 - t232) - t74 - s26 * (s26 * t241 + t17 * t149)) * s34 - 
     &s16 * (((-t188 - t73) * s25 - t200) * s25 - t201) - t50 * (t62 * t
     &210 + t238) - s26 * (s26 * t237 + t17 * t207)) * s34 - (((t231 - t
     &236) * s26 + t103 * t151) * s26 + t17 * t2) * s26 - t196) * s34 + 
     &t69 * (t34 * t47 + ddreal(71) * t12 + t215 * (((s16 * t214 + t153)
     & * s25 + ddreal(64) * t6) * s25 + ddreal(119) * t9) + t157 * t7 + 
     &t225 * t38 + t204 * t40 + t230 * t50 - ddreal(5) * t44 * t48)) * s
     &34
      t46 = t180 * (-t103 * (t77 + t222) + t17 * t40 + (s34 * t210 - t19
     &1) * s34 + t9) - t61 * (t6 * (((t55 + t3) * s25 + t81) * s25 + ddr
     &eal(14) * t9) + (s16 * (((t126 + t3) * s25 + t234) * s25 + t29) + 
     &((-t99 + (-t13 + t68) * s25 + t37) * s26 + ((t199 - t18) * s25 + t
     &81) * s25 + ddreal(62) * t9) * s26) * s26)
      t53 = t70 * t47
      t56 = t7 * t43
      t2 = t46 * s34 + t34 * s56 + t75 * (t6 * (((((-t135 - t115) * s25 
     &- t59) * s25 - t235) * s25 + ddreal(25) * t7) * s25 + t197) + ((((
     &-t192 + t166) * s26 + t151) * s26 + t2) * s26 + t196) * s26) + t76
     & * (-t103 * t49 + t6 * (((ddreal(47) * s16 + t3) * s25 + t113) * s
     &25 + ddreal(70) * t9) + (((s26 * t213 + t152) * s26 + t149) * s26 
     &+ t74) * s26) + t77 * (t6 * (((t38 * t5 - t156) * s25 - t212) * s2
     &5 - ddreal(70) * t7) - (s16 * (((t188 + t73) * s25 + t200) * s25 +
     & t201) + (((-s26 * t210 - t209) * s26 - t208) * s26 - t207) * s26)
     & * s26) + t79 * t22 + t80 * (t26 * t48 + t47 * (s26 * t226 + t224)
     & - t227 - t228 - t133 * s16 * ((t223 + t18) * s25 + ddreal(47) * t
     &6) - ddreal(52) * t165) + t60 + t141 * t6 * t182 + t53 * (t251 + t
     &184) * s16 - (-s26 + t135) * (s16 + t133) * t78 + t56
      t22 = -t6 + t38
      t20 = (t20 + t153) * s25 + t6
      t34 = ddreal(10) * s16
      t46 = t103 * t6
      t62 = (t34 + t153) * s25 + t46
      t44 = t17 * t44
      t73 = s16 * t103
      t74 = s26 * t129
      t81 = t26 * s16 * t25
      t83 = (t23 + t13) * t38
      t99 = t133 * t20
      t113 = (t25)**2
      t123 = t115 * t71
      t127 = s26 * s34
      t46 = ((-t44 * t76 - s34 * ((-s34 * t138 * t1 + t71) * s25 - t75 *
     & (t75 + t20))) * s56 - ((((-t103 * t25 + s34 - t13) * s34 + (t34 +
     & t74) * s25 + t54 * (s16 * s26 + t38) + t46) * s34 - t83 - t9 - t9
     &9) * s34 - t203 * ((s25 + t143 + t74) * s25 + t81)) * t75 - t123 *
     & t25 * s34) * s56 + (-t123 * t113 - t52 * t80) * s34 + (((((t133 +
     & t55 + t3) * s26 + t261) * s34 - t54 * (t248 + t222) - t17 * t62 *
     & s26 - t45) * s34 + t54 * t132 * t38 + (t17 * (t9 + t83) + t99) * 
     &s26 - ddreal(7) * t165) * s34 + t203 * t25 * ((t34 + t18 + t74) * 
     &s25 + t81)) * t75
      t74 = s16 + s26 - s34 + s56
      t59 = (t34 + s25) * s25 + t59
      t81 = ddreal(31) * s16
      t83 = ((ddreal(38) * s16 + s25) * s25 + t101) * s25 + t24
      t99 = ddreal(56) * s16
      t123 = t129 * t6
      t128 = (t99 + t198) * s25 + t123
      t139 = ddreal(78) * s16
      t144 = ddreal(23) * t7
      t143 = t143 + t187
      t146 = ddreal(92) * s16
      t147 = ddreal(102) * t9
      t148 = ((t146 + t98) * s25 + ddreal(91) * t6) * s25
      t149 = t147 + t148
      t151 = (t81 + t195) * s25 + ddreal(56) * t6
      t152 = ddreal(10) * t151
      t156 = ddreal(264)
      t146 = (((s25 * t156 - t146) * s25 - ddreal(504) * t6) * s25 - ddr
     &eal(484) * t9) * s25 + ddreal(912) * t7
      t156 = t214 * t9
      t157 = ((-t92 + t36) * s25 - t217) * s25 + t156
      t160 = ddreal(220) * t38 + ddreal(260) * t6
      t162 = ddreal(99)
      t164 = ddreal(84) * s25
      t166 = ddreal(252) * t9
      t167 = ((-ddreal(76) * s16 + t164) * s25 - ddreal(59) * t6) * s25
      t168 = t166 + t167
      t170 = ddreal(42) * s25
      t177 = (t170 + t94) * s25 + ddreal(72) * t6
      t180 = t13 + t94
      t182 = ddreal(10) * t177
      t184 = t168 * t5
      t185 = (-ddreal(40) * s25 + t211) * t38
      t186 = ddreal(320) * t6 * t22
      t187 = (t256 + s16) * s25 + t19
      t188 = t115 + t94
      t189 = -ddreal(55)
      t190 = s25 + t94
      t191 = t190 * t5
      t192 = ddreal(15) * t143 * s26
      t151 = ddreal(5) * t151
      t194 = ddreal(3) / ddreal(4) * t160 * s26
      t196 = ddreal(10) * t157
      t200 = ddreal(15) * t180 * s26
      t201 = ddreal(5) * t177
      t203 = t18 + t94
      t204 = ddreal(84) * t6
      t206 = ddreal(66) * s26
      t207 = t5 * t188
      t208 = t160 / ddreal(20)
      t84 = (((ddreal(130) * s16 + t84) * s25 + t202) * s25 + ddreal(235
     &) * t9) * s25 + t144
      t140 = ((((ddreal(102) * s16 + t13) * s25 + t140) * s25 + ddreal(2
     &11) * t9) * s25 + t169) * t103
      t169 = ((((s16 * t162 + t112) * s25 - ddreal(254) * t6) * s25 - dd
     &real(289) * t9) * s25 + ddreal(84) * t7) * s25 + ddreal(220) * t12
      t202 = t146 / ddreal(4)
      t209 = ((((-ddreal(80) * s16 + ddreal(63) * s25) * s25 - ddreal(23
     &4) * t6) * s25 - ddreal(357) * t9) * s25 + ddreal(389) * t7) * s25
     & + ddreal(75) * t12
      t210 = t169 * t103
      t211 = s25 * (t265 - t185) - t186
      t213 = ddreal(313)
      t216 = t211 * t103
      t219 = ddreal(301) * t6
      t213 = (s16 * t213 - t36) * s25
      t224 = ddreal(386) * t9
      t225 = ddreal(135) * t7
      t123 = ((s16 * t189 + t112) * s25 - t123) * s25 + ddreal(130) * t9
      t226 = -ddreal(203)
      t229 = ddreal(57) * s25
      t226 = ((s16 * t226 - t229) * s25 + ddreal(174) * t6) * s25 + ddre
     &al(145) * t9
      t230 = t123 * t54
      t126 = (s25 * t119 - t126) * s25 + ddreal(336) * t6
      t183 = (-t229 + t183) * s25 + t205 * t6
      t205 = ddreal(33) * s16
      t229 = s16 * t41
      t231 = t9 * t40
      t232 = t7 * s25
      t234 = s16 * (t205 - t91)
      t235 = ((-t219 - t213) * s25 + t224) * s25 + t225
      t190 = (t141 * t83 + t229) * s16 + t79 * ((t243 * t190 + t234) * s
     &26 + ddreal(8) * t45)
      t34 = t17 * t190 - ((((-t129 * t187 * t50 - (-t47 * t126 - t215 * 
     &t183 - ddreal(24) * t203 * t50 - ddreal(98) * t165 + t178) * s34 +
     & t102 * t103 - ddreal(20) * t188 * t48 + t204 * t40 - t215 * t226 
     &- t230 * t47 - ddreal(250) * t232) * s34 - t155 * t125 + t173 * t1
     &80 + t184 * t50 + t201 * t48 + t215 * t235 + t216 * t47 + ddreal(3
     &40) * t171 - ddreal(188) * t231 - ddreal(14) * t30) * s34 - t210 *
     & t47 - t196 * t48 - t90 * t43 + ddreal(84) * t12 * t38 + ddreal(12
     &8) * t9 * t39 + ddreal(138) * t7 * t40 - t150 * t11 - t215 * t209 
     &- t50 * (ddreal(3) / ddreal(10) * t160 * t47 + t146)) * s34 + t69 
     &* (t49 * t54 * t143 + t5 * t149 * t50 + t231 * t125 + t140 * t47 +
     & t155 * t145 + t151 * t48 + t215 * t84 + ddreal(106) * t171 + ddre
     &al(104) * t30)) * s34 + ddreal(18) * t174 + ddreal(18) * t175 - t3
     &4 * s26 * t78 + t131 * t59 * t50 + t254 * t137 * t48 + ddreal(106)
     & * t12 * t39 - t56 * t51 - ddreal(41) * t10 * t38 + ddreal(23) * t
     &60 - ddreal(14) * t11 * t40 + t163 * s16 * t128 * t47
      t36 = -(s16 * t235 + (((t180 * s26 + t177) * s26 + t168) * s26 + t
     &211) * s26) * t47 + t45 * (s16 * (((t195 + t233) * s25 - t51 * t6)
     & * s25 - ddreal(40) * t9) - t52 * (((-s16 * t51 - t36) * s25 - t42
     &) * s25 + t240))
      t19 = (t103 * t89 - t129 * t88 + t228 - ((t180 * s34 - t208) * s34
     & + t143 * t69) * s34) * s56 + s25 * ((t40 * t54 + t156) * s25 - t2
     &27) - t129 * s16 * t181 - (((s16 * (t220 + t66) + (ddreal(24) * s2
     &6 + t94 + t170) * s25) * s34 - ((s25 * t125 + t206 - t233) * s25 -
     & t19) * s25 - t6 * (ddreal(96) * s16 + ddreal(78) * s26)) * s34 + 
     &t69 * (s16 * (t99 + t206) + (t120 + t81 + t195) * s25)) * s34 + t2
     &07 * t76 + ddreal(18) * t141
      t19 = t19 * s56 + s25 * ((((-t204 + t106) * s25 + ddreal(96) * t9)
     & * s25 - t26 * t7) * s25 - t197) - (((-(s16 * (ddreal(138) * s16 +
     & ddreal(100) * s26) + (s25 * t193 + t120 + t73) * s25) * s34 + (t2
     &00 + t201) * s26 + t166 + t167) * s34 - (((-t35 + t104) * s25 - dd
     &real(126) * t6) * s25 - ddreal(121) * t9) * s25 - (t196 + t194) * 
     &s26 - ddreal(228) * t7) * s34 + t69 * ((t151 + t192) * s26 + t147 
     &+ t148)) * s34 + ddreal(18) * s16 * (t12 + t43) + ddreal(30) * t14
     &1 * t137 - t161 * t203 + ddreal(45) * t53
      t19 = t19 * s56 + t119 * t53 * (s26 + t137) + t129 * (t141 * t59 +
     & t10) + t76 * (t17 * (s16 * ((s25 * t189 - t159) * s25 + ddreal(13
     &0) * t6) + ddreal(20) * t188 * t47 + ddreal(20) * t40 + t243 * t18
     &7) + t191 * t75) - (((t75 * (s16 * (ddreal(112) * s16 + ddreal(120
     &) * s26) + (s26 * t214 + t112 - t85) * s25) - s25 * (-t265 + t185)
     & - ((-ddreal(20) * t180 * s26 - t182) * s26 - t184) * s26 - t186) 
     &* s34 - s16 * ((((s25 * t162 - t136) * s25 - ddreal(289) * t6) * s
     &25 + ddreal(84) * t9) * s25 + ddreal(220) * t7) - ddreal(20) * t15
     &7 * t47 - ddreal(20) * t43 - s26 * (t160 * t47 + t146)) * s34 + t6
     &9 * (s16 * (((ddreal(110) * s16 + t124) * s25 + ddreal(211) * t6) 
     &* s25 + t212) + t47 * (ddreal(20) * t143 * s26 + t152) + t5 * (s26
     & * t149 + t39))) * s34 + t87 * (((((-ddreal(20) * s16 + t198) * s2
     &5 - ddreal(98) * t6) * s25 + ddreal(212) * t9) * s25 - ddreal(133)
     & * t7) * s25 + ddreal(8) * t12)
      t10 = (t19 * s56 + s16 * (-ddreal(5) * t78 + t87 * (((((s25 * t72 
     &- t32) * s25 + ddreal(35) * t6) * s25 + ddreal(125) * t9) * s25 - 
     &ddreal(136) * t7) * s25 + t12 * t145)) + t103 * s16 * (t141 * t128
     & + t10) - (((((-((t120 + t205 - t91) * s16 + t142 * s26) * s34 + s
     &16 * t183 + (t176 * t203 + t126) * s26) * s34 - s16 * t226 - ((t26
     &4 * t188 + ddreal(18) * t187) * s26 + t230) * s26) * s34 - s16 * (
     &((t219 + t213) * s25 - t224) * s25 - t225) - (((-t182 - t200) * s2
     &6 - t168 * t54) * s26 - t216) * s26) * s34 - s16 * t209 - (((ddrea
     &l(20) * t157 + t194) * s26 + t202 * t54) * s26 + t210) * s26) * s3
     &4 + t69 * (s16 * t84 + (((t152 + t192) * s26 + t149 * t54) * s26 +
     & t140) * s26)) * s34 + t229 + ddreal(18) * t53 * t59 + t119 * t70 
     &* t137 * t50 + ddreal(45) * t70 * t48) * s56 + t34
      t12 = t61 * (((((t143 * s26 + (t199 + t195) * s25 + t154) * s26 + 
     &((ddreal(62) * s16 + t98) * s25 + ddreal(70) * t6) * s25 + t29) * 
     &s26 + (((t139 + t13) * s25 + t101) * s25 + t31) * s25 + t144) * s2
     &6 + t87 * (((t218 + t249) * s25 + t221) * s25 + t130)) * s26 + t45
     & * (((t135 + t21) * s25 + t15) * s25 + t111)) * s34
      t8 = t10 * s56 + t38 * t8 + t41 * t9 + t75 * ((s16 * t209 + (((s26
     & * t208 + t17 * t157) * s26 + t202) * s26 + t169) * s26) * t47 + t
     &45 * (t6 * (((-t139 - t95) * s25 + ddreal(66) * t6) * s25 + t265) 
     &+ t52 * ((((t218 - t247) * s25 - t172) * s25 - t108) * s25 + ddrea
     &l(130) * t7))) + t76 * (-t103 * (-t187 * t48 + t231) + s26 * ((t17
     & * t123 * s26 + s16 * t226) * s26 + t45 * ((-t164 - t73) * s25 + d
     &dreal(250) * t6)) + t207 * t49 + t232 * (ddreal(61) * s25 + t223))
     & + t77 * t36 + t79 * (((t191 * s26 + t234) * s26 + ddreal(16) * t4
     &5) * s26 + t165 * t17) + t80 * (-t52 * t45 * (s16 * t96 - t115) - 
     &t54 * t203 * t48 - ddreal(15) * t102 - t142 * t7 - t47 * (s16 * t1
     &83 + s26 * t126 / ddreal(3))) - t12 + t53 * (s26 * t128 + t83) * s
     &16 + t141 * t45 * ((t81 + t18) * s25 + t37)
      t10 = ddreal(1) / s16
      t12 = ddreal(1) / (t64)**2
      t15 = ddreal(1) / (t63)**2
      t19 = ddreal(1) / (t58)**2
      t2 = t27 * MB1001(0) + ((t103 * t70 * t48 * (t47 + t59) - ddreal(5
     &) * s16 * (t47 * t78 - t60) + t8 + t173 * t70 * t137 + t171 * (t38
     & * ((t242 - t250) * s25 - t217) + t97)) * MB1110(1) - t17 * (t70 *
     & t4 * t48 + t41 * t6) - t5 * t11 * (t9 + t40) - ddreal(7) * t174 +
     & ddreal(7) * t175 - t14 + t134 * t47 - t30 * (t38 * (t23 - t21) - 
     &t24) - t28 * (s16 * (t7 * ((s25 * t65 - t24) * s25 + t97) + t114) 
     &+ t2) * MB1011(1)) * t10
      t4 = (t28)**2
      t1 = t10 * (t12 * ((s34 * ((-t71 * t25 * t113 + s34 * t113 * ((s26
     & * t5 + s25) * s25 + t103 * (s25 + s26 + s16) * s16) * t1) * s25 +
     & t75 * (s34 * (t87 * t100 * s26 - t44 * t50 - t47 * t62 - t118 + t
     &165 + t127 * ((t73 + s26 + t13) * s26 + t261 - t127)) - t25 * (t17
     & * t222 * t1 + (-s26 * t20 - t13 * (t261 + t22)) * s26 + t165 * t1
     &03))) + t46 * s56) * MB1101(0) * t19 - t28 * t4 * t74 * MB1111D2(0
     &)) + (t74)**2 * MB1110(0) * t15)
      result = t12 * t19 * (-(-t103 * (t47 * (t80 - t141) - (-t5 * t88 +
     & (-t138 * t69 + t75 * (-t117 + s34)) * s34 + t39 + t7 + t178 + t17
     & * t107 * t75) * s56 * t179) - t17 * t9 * t181 + t82 * s56 + t75 *
     & (t6 * ((ddreal(5) * t6 + t110) * s25 + t111) + s26 * (s26 * t105 
     &+ t109) + t54 * t107 * t50) + t76 * (((t133 + t122 + t55) * s26 + 
     &t132 * t54) * s26 + t17 * t9) - t77 * (t103 * (t117 * t50 + t118) 
     &+ t9 * (t153 + t68) + s26 * (s26 * t116 + t121)) - t61 * (ddreal(8
     &) * t6 * t16 + s26 * (t93 * s26 + (t92 + t91) * s25 + t57) + t94 *
     & t38) * s34 + t33 * (t18 + t67) + t70 * t86 * t47 + t158 + t171 * 
     &(-t67 + t3)) * I300s16(1) + t15 * t2) / ddreal(12) + t1 / ddreal(6
     &)

           intHs16s25s26s34s56x1411D6eps0 = result
       end function intHs16s25s26s34s56x1411D6eps0

       function intHs16s25s26s34s56x1411D6eps1()
           implicit none
           type(dd_complex) :: intHs16s25s26s34s56x1411D6eps1
           type(dd_complex) ::  t1,t2

           type(dd_complex) :: result

      t1 = s16 + s26 - s34 + s56
      t2 = s16 + s25 + s26 - s34 + s56
      t2 = ddreal(1) / (t2)**2
      result = ddreal(1) / s16 * t2 * (-s25 / ddreal(12) - t1 / ddreal(4
     &) + (t1)**2 * MB1110(1) / ddreal(6))

           intHs16s25s26s34s56x1411D6eps1 = result
       end function intHs16s25s26s34s56x1411D6eps1


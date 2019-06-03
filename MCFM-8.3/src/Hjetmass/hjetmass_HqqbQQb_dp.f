      double complex function hjetmass_qqbqqb_bubble_pmpm_mhsq_dp
     &  (i1,i2,i3,i4,za,zb)
      implicit double complex (t)
      integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      double complex ret
      double complex hjetmass_qqbqqb_bubble_pmpm_s12_dp
      double complex hjetmass_qqbqqb_bubble_pmpm_s34_dp

      hjetmass_qqbqqb_bubble_pmpm_mhsq_dp =
     &   - hjetmass_qqbqqb_bubble_pmpm_s12_dp(i1,i2,i3,i4,za,zb)
     &   - hjetmass_qqbqqb_bubble_pmpm_s34_dp(i1,i2,i3,i4,za,zb)
      return
      
      end function

      double complex function hjetmass_qqbqqb_bubble_pmpm_s12_dp
     &  (i1,i2,i3,i4,za,zb)
      implicit double complex (t)
      integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      double complex ret

      t1 = za(i1, i3)
      t2 = za(i1, i4)
      t3 = za(i2, i3)
      t4 = za(i2, i4)
      t5 = zb(i3, i1)
      t6 = zb(i3, i2)
      t7 = zb(i4, i2)
      t8 = zb(i4, i1)
      t9 = za(i3, i4)
      t10 = zb(i4, i3)
      t11 = t5 ** 2
      t12 = t11 ** 2
      t13 = t5 * t12
      t14 = t5 * t11
      t15 = t4 ** 2
      t16 = t15 ** 2
      t17 = t4 * t16
      t18 = t4 * t15
      t19 = t3 ** 2
      t20 = t19 ** 2
      t21 = t3 * t19
      t22 = t8 ** 2
      t23 = t22 ** 2
      t24 = t8 * t22
      t25 = t15 * t22
      t26 = t19 * t11
      t27 = t26 + t25
      t28 = t1 ** 2
      t29 = t28 ** 2
      t30 = t1 * t29
      t31 = t1 * t28
      t32 = t2 ** 2
      t33 = t32 ** 2
      t34 = t2 * t33
      t35 = t2 * t32
      t36 = t6 ** 2
      t37 = t36 ** 2
      t38 = t6 * t37
      t39 = t6 * t36
      t40 = t1 * t4
      t41 = t40 * t2
      t42 = t3 * t5
      t43 = t42 * t32
      t44 = t28 * t5
      t45 = t4 * t8
      t46 = t45 * t28
      t47 = t20 * t11
      t48 = t47 * t36
      t49 = t1 * t5
      t50 = t2 * t24
      t51 = t19 * t14
      t52 = t51 * t1
      t53 = t44 * t6
      t54 = t32 * t23
      t55 = t54 * t15
      t56 = t28 * t36
      t57 = t51 * t7
      t58 = t15 * t6
      t59 = t3 * t6
      t60 = t59 * t18
      t61 = t1 * t32
      t62 = t50 * t28
      t63 = t1 * t6
      t64 = t2 * t7
      t65 = t64 * t63
      t66 = t5 * t4
      t67 = t7 ** 2
      t68 = t67 ** 2
      t69 = t7 * t67
      t70 = t28 * t18 * t7
      t71 = t1 * t16
      t72 = t17 * t69
      t73 = t31 * t15
      t74 = t73 * t8
      t75 = t21 * t39
      t76 = t75 * t18
      t77 = t35 * t24
      t78 = t77 * t51
      t79 = t63 * t2
      t80 = t64 * t17
      t81 = t28 * t11
      t82 = t81 * t58
      t83 = t63 * t7
      t41 = t83 * (t2 * (t3 * (t71 * t67 * t8 * t14 - t72 * t8 * t11 + t
     #70 * t8 * t12 - t74 * t13) + t76 * t5 * t22 + t78 * t4 - t79 * t20
     # * t13 - t80 * t6 * t23) + t82 * (t26 + t25) * t10 * t9) + t66 * (
     #t7 * (t7 * (t3 * (-t44 * t18 * t36 * t8 - t46 * t32 * t14) + t7 * 
     #(t41 * (t5 * t27 * t2 - t4 * t27 * t6) - t43 * t18 * t8 * t7) + t4
     #1 * t21 * t11 * t36 + t33 * t15 * t23 + t48 * t32) + t53 * t2 * (t
     #15 * (-t49 * t22 + t50) - t52)) + t56 * (-t46 * t3 * t14 + t48 + t
     #55)) * t10 * t9 + t65 * (-t62 * t18 * t14 + t60 * t7 * (t58 * t24 
     #- t57) + t61 * t21 * t13 * t8)
      t46 = t64 - t63
      t48 = t2 * t8
      t84 = t48 - t59
      t85 = t48 + t59
      t86 = -t64 + t63
      t87 = t15 * t67
      t88 = t87 * t59
      t89 = t18 * t6
      t90 = t89 * t7 * t8
      t91 = t3 * t2
      t92 = t18 * t69
      t93 = t18 * t24
      t94 = t32 * t22
      t95 = t32 * t21
      t96 = t64 * t19
      t97 = t32 * t67
      t98 = t56 + t97
      t99 = t29 * t11
      t100 = t35 * t69
      t101 = t31 * t39
      t102 = t19 * t36
      t103 = t91 * t1
      t104 = t103 * t14
      t105 = t64 * t1
      t106 = t31 * t14
      t107 = t32 * t24
      t46 = t4 * (t18 * (-t106 * t6 * t8 * t85 * t67 + t104 * (-t48 + t5
     #9) * t68 - t105 * t36 * t24 * (t102 + t81) + t77 * t59 * t5 * t69)
     # + t91 * t31 * t13 * t84 * t67 * t4 - t107 * t63 * t16 * t5 * t69 
     #+ t79 * t17 * t5 * t22 * t68 + t5 * (t26 * (t81 * t98 + t87 * t98)
     # + t14 * (-t101 + t100) * t21 + t93 * (t101 - t100) + t25 * (t32 *
     # t15 * t68 + t99 * t36)) * t10 * t9 - t78 * t56) + t5 * (t1 * (t87
     # * t47 * t2 * t39 - t95 * t6 * t7 * (t94 + t87) * t14) + t28 * (-t
     #92 * t91 * t85 * t14 + t90 * (-t48 * t87 + t77 + t88) * t5) + t29 
     #* (t96 * t4 * t6 * t13 + t90 * t84 * t14) + t31 * (-t64 * t21 * t4
     # * t36 * t12 + t75 * t45 * t2 * t14) + t65 * t4 * (t46 * t14 * t21
     # + t25 * t42 * t86 + t26 * t45 * t46 + t93 * t86) * t10 * t9 - t94
     # * t76 * t67)
      t76 = t56 + t97
      t78 = t1 * t11
      t79 = t105 * t59
      t84 = t79 * t15 * t8
      t86 = t15 * t8
      t90 = t86 * t35 * t68
      t105 = t18 * t22
      t108 = t32 * t19
      t109 = t18 * t67
      t110 = t28 * t3
      t111 = t97 * t3
      t112 = t1 * t2
      t113 = t4 * t7
      t114 = t36 * t15
      t115 = t11 * t3
      t46 = t8 * (t22 * (t113 * (t16 * (-t110 * t37 - t111 * t36) - t112
     # * t17 * t36 * t67 - t57 * t34) - t71 * t35 * t36 * t7 * t22) + t1
     #15 * (t32 * t17 * t7 * t68 - t51 * t31 * t32 * t36 + t114 * t30 * 
     #t14)) + t5 * (t5 * (t5 * (t18 * (t6 * t67 * t76 * t21 + t62 * t76)
     # + t78 * t21 * (-t91 * t28 * t39 - t33 * t67 * t8 + t101 * t4)) - 
     #t84 * (t75 + t77)) + t105 * (-t28 * t21 * t38 + t90)) + t46 + t5 *
     # t15 * (t5 * (t5 * (t7 * (t7 * (t108 * (t48 - t59) * t7 + t94 * t2
     #8 * t4) - t75 * t28) - t101 * t19 * t8) + t22 * t4 * (t3 * (t101 -
     # t100) - t48 * t1 * t76)) + t109 * t56 * t22) * t10 * t9 - t72 * t
     #35 * t6 * t23 - t95 * t83 * t12 * (t102 + t81)
      t62 = t1 * t21
      t71 = t20 * t12
      t72 = t16 * t6
      t75 = t2 * t21
      t77 = t113 * t59
      t95 = t48 * t49
      t101 = t6 * t2
      t116 = t115 * t86
      t117 = t26 * t58 * t2
      t118 = t7 * t5
      t119 = -t118 * t19 + t40 * t22
      t120 = t115 * t87
      t121 = t66 * t1
      t122 = t121 * t7
      t119 = t92 * t11 * t8 * t119 * t35 - t73 * t3 * t14 * t119 * t39 -
     # t33 * t20 * t13 * t69 - t31 * t17 * t37 * t23 + t122 * (t6 * (t6 
     #* (t86 * t81 * t3 + t101 * t47) + t4 * t2 * (t8 * (t8 * (t49 * t15
     # * t7 + t94 * t4) - t42 * (t81 + t87)) + t57 * t1)) + t120 * t32 *
     # t8) * t10 * t9
      t114 = t11 * t32 + t114
      t82 = t82 * t2 * t67
      t123 = t66 * t10 * t9 * t8 * t3
      t124 = t26 * t7
      t125 = t66 * t7
      t50 = t50 * t15
      t126 = t11 * t6 * t21
      t127 = t48 * t66
      t54 = t15 * (t21 * t11 * t36 * t9 * t10 * t76 * t8 + t115 * t107 *
     # t9 * t10 * t76) + t54 * t60 * t5 * t76 + t118 * t59 * t9 * t10 * 
     #t76 * t22 * t16 + t51 * t45 * t2 * (t49 * t76 * t10 * t9 + t102 * 
     #t76) - t80 * t28 * t39 * t23 - t63 * t35 * t20 * t13 * t67
      t80 = t21 * t35
      t128 = t80 * t12
      t129 = t1 * t36
      t130 = t8 * t3
      t120 = t130 * (t7 * (t120 * t35 * (-t45 * t6 + t78) + t129 * (t128
     # + t72 * (t81 * t8 + t107))) - t91 * t74 * t14 * t39 + t65 * t15 *
     # t11 * (t95 + t77) * t10 * t9)
      t131 = t1 * t7
      t132 = t91 * t11
      t51 = t130 * t66 * (t73 * t2 * t5 * t39 * t22 + t100 * t26 * t4 * 
     #t6 + t105 * t31 * t37 + t51 * t33 * t69)
      t73 = t131 * t130 * t18 * t14 * (t67 * t8 * t35 + t110 * t39)
      t70 = t66 * (t8 * (t8 * (t3 * (t2 * (t2 * (t2 * (-t28 * t4 * t14 *
     # t67 + t104 * t67) + t3 * t36 * (t106 + t92)) - t29 * t4 * t14 * t
     #36) + t70 * t3 * t37) - t29 * t18 * t5 * t39 * t8) - t89 * t19 * t
     #5 * t67 * t98) - t80 * t4 * t14 * t68)
      t89 = t94 + t87
      t100 = -t40 + t91
      t104 = t67 * t22
      t106 = t107 * t16 * t67
      t107 = t112 * t21
      t112 = t15 * t7
      t15 = t5 * (t5 * (t5 * (t5 * (t5 * (-t108 * t66 * t31 * t67 + t112
     # * t28 * t19 * t76) + t61 * (t18 * (-t104 * t28 + t19 * t68) + t10
     #7 * t36 * t22)) + t112 * (-t87 * t56 * t19 - t36 * t76 * t20 + t25
     # * t28 * t76)) + t86 * (t109 * t1 * t8 * t76 + t34 * t3 * t67 * t2
     #4 - t61 * t24 * t76 * t4 + t28 * t20 * t38)) - t56 * t15 * t16 * t
     #69 * t22) + t6 * (t6 * (t11 * (t110 * t35 * t15 * t23 + t26 * (t18
     # * t31 + t80) * t67 - t105 * t30 * t11 + t29 * t19 * t12 * t100) +
     # t6 * (t6 * (t72 * t28 * t19 * t24 + t71 * t28 * t100) + t86 * (t2
     #8 * (t104 * t16 + t55) + t97 * t47))) + t106 * t89) + t97 * t5 * (
     #t21 * (t28 * t2 * t13 + t35 * t22 * t14) - t22 * t16 * t7 * t89 - 
     #t26 * t16 * t69) + t106 * t19 * t39
      t29 = t125 * t61 * t36
      t30 = t48 + t59
      t34 = t56 * t30 + t97 * t30
      t34 = t122 * (t36 * (t127 * t19 * (t81 - t87) + t25 * t91 * t81) +
     # t39 * (t28 * t16 * t24 + t105 * t96) + t6 * (t43 * (t103 * t11 - 
     #t81 * t4 + t109) * t22 + t26 * t87 * t32 * t8) + t128 * t67) + (t7
     # * t6 * t24 * t76 * t17 + t107 * t76 * t13 + t75 * t34 * t12 + t72
     # * t24 * t34) * t10 * t9
      t11 = t84 * t11 * (t102 * t49 + t113 * t94)
      t38 = t29 * t130
      t43 = t25 * t26
      t47 = t43 * t65 * (t118 * t32 + t40 * t36)
      t55 = t124 * t25 * t61 * t36
      t26 = t38 * (t8 * (-t26 * t2 + t86 * (t49 - t59)) + t113 * t26)
      t57 = t79 * t66 * t8 * (t129 * t105 + t57 * t32)
      t68 = t48 * t88 * t81 * (t58 * t8 + t132)
      t79 = t55 * (t113 + t49)
      t80 = t116 * (t99 * t3 * t39 + t90)
      t43 = t43 * (t31 * t4 * t37 + t5 * t33 * t69)
      t11 = -88 * t55 * t9 * t10 - t15 + 64 * t26 + 69 * t79 + 78 * t47 
     #+ 30 * t38 * (t105 * t7 + t52) + 50 * t11 + 54 * t57 + 7 * t80 + 2
     #1 * t29 * (t1 * t18 * t23 + t14 * t7 * t20) + 26 * t51 + 29 * t43 
     #- 3 * t70 - 6 * t34 - 13 * t82 * (t126 + t50) - 22 * t73 - 33 * t6
     #8 - 32 * t101 * t66 * t22 * t19 * (t131 * t114 * t10 * t9 + t58 * 
     #t98 * t8 + t132 * t98) - 28 * t127 * t59 * ((t109 * t63 * t8 + t1 
     #* (t19 * (t125 + t78) + t25 * t1) * t36 + t7 * (t121 * t22 + t25 *
     # t7 + t124) * t32 + t110 * t64 * t14) * t10 * t9 + t65 * (t126 + t
     #50)) - 36 * t48 * t60 * t28 * t14 * t67 * t30 - 35 * t117 * t22 * 
     #(t56 * t85 + t97 * t85) - 70 * t55 * t30
      t15 = -12
      t18 = -10
      t19 = t59 - t48 + t113 - t49
      t20 = -2 * t48 * t59 + 2 * t49 * (t59 + t48 - t113) + 2 * t113 * t
     #30
      t25 = 4 * t42 * t64 + 4 * t45 * t63
      t26 = t20 + t25 + t102 + t81 + t87 + t94
      t28 = t19 ** 2
      t29 = (4 * t45 + 4 * t42) * (t63 + t64)
      t19 = t29 - t19 * cdsqrt(t26) + t28
      t26 = t29 + t28
      t26 = cdsqrt(t26)
      t20 = -t20 - t25 + t48 * t26 + t49 * t26 - t59 * t26 - t102 - t81 
     #- t94 - t113 * (t113 + t26)
      t25 = 0.1D1 / t10
      t26 = 0.1D1 / t9
      t28 = 256
      t19 = 0.1D1 / t19 ** 2
      t20 = 0.1D1 / t20 ** 2
      ret = t28 * (-4 * t41 + 2 * t46 + t15 * (t2 * (t1 * t17 * t36 * t2
     #4 * t9 * t10 * t67 + t40 * t3 * t39 * (t93 * t9 * t10 - t62 * t12)
     # * t7) + t35 * (t83 * t21 * t12 * t8 * t9 * t10 - t63 * t16 * t5 *
     # t23 * t67) - t74 * t21 * t14 * t37 - t115 * t92 * t33 * t24 + t83
     # * t10 * t9 * (t62 * t13 + t72 * t23 + t71 * t6) * t32) + t18 * (t
     #116 * (t77 * t76 + t95 * t76) * t10 * t9 + t101 * (t8 * (t22 * (t5
     #3 * t17 * t67 + t72 * (-t31 * t5 * t6 + t56 * t3 + t111) * t8) + t
     #71 * t2 * t76) + t75 * t4 * t12 * t67 * (-t64 * t3 + t44))) - 44 *
     # t117 * t9 * t22 * t10 * t76 - 8 * t119 - 16 * t123 * (t82 + (t56 
     #* t114 + t97 * t114) * t8 * t3) - 14 * t54 - 56 * t123 * t61 * t36
     # * t7 * t27 - 20 * t120 + t11) * t26 * t19 * t20 * t25

      hjetmass_qqbqqb_bubble_pmpm_s12_dp =  ret/16d0*(0d0,1d0)
      return
      
      end function

      double complex function hjetmass_qqbqqb_bubble_pmpm_s34_dp
     &  (i1,i2,i3,i4,za,zb)
      implicit double complex (t)
      integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      double complex ret

      t1 = za(i1, i3)
      t2 = za(i1, i4)
      t3 = za(i2, i3)
      t4 = za(i2, i4)
      t5 = zb(i3, i1)
      t6 = zb(i3, i2)
      t7 = zb(i4, i2)
      t8 = za(i1, i2)
      t9 = zb(i2, i1)
      t10 = zb(i4, i1)
      t11 = t10 ** 2
      t12 = t11 ** 2
      t13 = t10 * t11
      t14 = t1 ** 2
      t15 = t14 ** 2
      t16 = t1 * t14
      t17 = t3 ** 2
      t18 = t17 ** 2
      t19 = t3 * t17
      t20 = t7 ** 2
      t21 = t20 ** 2
      t22 = t7 * t20
      t23 = t14 * t11
      t24 = t17 * t20
      t25 = t24 + t23
      t26 = t5 ** 2
      t27 = t26 ** 2
      t28 = t5 * t27
      t29 = t5 * t26
      t30 = t2 ** 2
      t31 = t30 ** 2
      t32 = t2 * t30
      t33 = t4 ** 2
      t34 = t33 ** 2
      t35 = t4 * t34
      t36 = t4 * t33
      t37 = t2 * t3
      t38 = t37 * t1
      t39 = t14 * t4
      t40 = t36 * t22
      t41 = t2 * t11
      t42 = t14 * t2
      t43 = t42 * t36
      t44 = t15 * t36
      t45 = t30 * t10
      t46 = t36 * t20
      t47 = t32 * t29
      t48 = t5 * t4
      t49 = t3 * t6
      t50 = t2 * t10
      t51 = -t49 + t50
      t52 = t6 ** 2
      t53 = t52 ** 2
      t54 = t6 * t52
      t55 = t30 * t26
      t56 = t33 * t52
      t57 = t56 + t55
      t58 = t36 * t54
      t59 = -t58 + t47
      t60 = t32 * t18 * t6 * t7
      t61 = t32 * t3
      t62 = t17 * t22
      t63 = t62 * t1 * t2
      t64 = t10 * t7
      t65 = t30 * t11
      t66 = t4 * t7
      t67 = t3 * t7
      t68 = t67 * t45
      t69 = t6 * t1
      t70 = t61 * t4 * t20
      t71 = t1 * t2
      t72 = t71 * t19
      t73 = t46 * t14 * t6
      t74 = t1 * t4
      t75 = t17 * t6
      t76 = t17 * t52
      t77 = t19 * t22
      t78 = t16 * t13
      t79 = t30 * t12
      t80 = t1 * t5
      t81 = t2 * t17
      t82 = t81 * t54
      t83 = t74 * t52
      t84 = t14 * t30
      t85 = t47 * t1
      t86 = t1 * t33
      t87 = t5 * t20
      t88 = t52 * t11
      t89 = t19 * t20
      t90 = t30 * t17
      t91 = t90 * t5
      t92 = t14 * t27
      t93 = t48 * (t66 * (t10 * (t3 * (t7 * (t7 * (t4 * (t4 * (t83 * (t4
     #9 + t80) - t82) + t84 * t29) + t85 * t3 - t86 * (t56 + t55) * t7) 
     #+ t56 * t16 * t26) + t72 * t5 * t53) - t87 * t43 * t6 * t10) + t63
     # * t33 * t26 * t6) + (t4 * (t4 * (t5 * (t1 * (t1 * (t88 * (-t50 + 
     #t49) * t1 - t87 * (t76 + t65)) + t89 * t54) - t91 * t21) - t56 * t
     #20 * t25) + t55 * t6 * (t78 - t77)) - t90 * t92 * t20) * t9 * t8)
      t94 = t14 * t26
      t95 = t33 * t20
      t96 = t95 + t94
      t97 = t14 * t5
      t98 = t94 * t37
      t99 = t1 * t19
      t100 = t3 * t35
      t101 = t7 * t17
      t102 = t14 * t13
      t103 = t2 * t5
      t104 = t103 * t34
      t105 = t24 * t48 * t6 * (t55 * t19 * t52 - t56 * t96 * t3 - t104 *
     # t22)
      t106 = t33 * t7
      t107 = t106 * t26
      t108 = t107 * t2 * ((t24 * t10 + t102) * t5 * t30 + t89 * t4 * t52
     #) * t9 * t8
      t109 = t1 * t26
      t110 = t64 * t2
      t111 = t1 * t7
      t112 = t111 * (t17 * (-t66 * t32 * t8 * t9 * t10 * t27 - t106 * t4
     #5 * t6 * t8 * t9 * t29 + t110 * t36 * t52 * t8 * t9 * t26 + t64 * 
     #t34 * t54 * t8 * t9 * t5) + t3 * (-t109 * t41 * t36 * t52 * t8 * t
     #9 + t74 * t32 * t8 * t9 * t11 * t27 - t80 * t34 * t54 * t8 * t9 * 
     #t11 + t65 * t69 * t33 * t8 * t9 * t29) + t60 * t28 + t71 * t35 * t
     #54 * t12)
      t51 = t112 + t10 * (t10 * (t10 * (t10 * (t97 * t32 * t36 * t52 * t
     #10 + t67 * t69 * t31 * t33 * t26) + t85 * (t91 * t7 - t74 * t96 + 
     #t98)) + t6 * (-t44 * t3 * t29 * t52 + t99 * t34 * t53 * t7 + t100 
     #* t52 * (t74 + t37) * t22 + t16 * t2 * (-t14 * t33 + t90) * t28)) 
     #+ t101 * (t32 * (t16 * t26 * t27 - t40 * t29) + t100 * t53 * t20))
     # + t93 + t5 * (t5 * (t5 * (t5 * (t5 * (-t68 * t48 * t15 + t66 * t3
     #7 * t16 * (t66 * t51 + t65)) + t69 * (t64 * t51 * t36 * t16 - t61 
     #* t39 * t13 + t60 * t10 + t63 * t36)) + t74 * t41 * (-t70 * t10 + 
     #t72 * t54 + t73)) + t58 * t64 * t14 * t3 * (t74 * t10 - t75)) + t7
     #6 * t46 * t32 * t13) + t48 * (-t26 * t57 * t11 * t15 - t77 * t59 +
     # t78 * t59) * t9 * t8 + t79 * t67 * t1 * t34 * t54 + t105 + t108
      t59 = t24 + t23
      t60 = t50 * t59
      t63 = t49 * t59 + t60
      t72 = t38 * t26
      t85 = t81 * t5
      t91 = t89 * t32
      t93 = t10 * t54
      t100 = t111 * t48
      t61 = t100 * (t52 * (t67 * t30 * t36 * t13 + t85 * (-t39 * t26 + t
     #46 + t72) * t10 + t94 * t41 * t3 * t33) + t6 * (t24 * t45 * t33 * 
     #t26 + t65 * t48 * t3 * (-t95 + t94)) + t91 * t27 + t102 * t34 * t5
     #4) + (t61 * t1 * t59 * t28 + t64 * t54 * t59 * t35 + t61 * t63 * t
     #27 + t93 * t63 * t34) * t9 * t8
      t63 = t33 * t11
      t105 = t17 * t26
      t108 = t37 * t26
      t112 = t45 * t48 * t52 * t3 * (t108 * t25 + t111 * (t105 + t63) * 
     #t9 * t8 + t10 * t25 * t6 * t33)
      t113 = t37 * t29
      t99 = t99 * t32
      t114 = t69 * t33
      t115 = t6 * t2
      t116 = t69 * t3 * t5
      t107 = t115 * (t10 * (t10 * (-t99 * t27 * t7 + t114 * (t14 * (t113
     # - t107) - t106 * t76) * t10) + t115 * t77 * t36 * t26) - t86 * t6
     #4 * t3 * t26 * (t50 * t66 + t116) * t9 * t8 - t71 * t77 * t33 * t2
     #7)
      t117 = t56 * t45 * t26 * t3
      t118 = t115 * t48
      t119 = t118 * (t36 * t16 * t52 * t12 + t30 * t18 * t29 * t22 + t77
     # * t45 * t4 * t26 + t56 * t78 * t3 * t5)
      t120 = t76 + t95
      t121 = t14 * t3
      t122 = t16 * t4
      t123 = t90 * t34
      t124 = t32 * t10
      t125 = t11 * t59
      t126 = t125 * t31
      t127 = t69 * t46
      t128 = t33 * t6
      t75 = t5 * (t5 * (t5 * (t5 * (t1 * (t36 * (t90 * t21 + t14 * (-t76
     # + t65) * t20) + t99 * t88) + t5 * (-t90 * t48 * t16 * t20 + t65 *
     # t106 * t15)) + t106 * (-t95 * t84 * t11 + t56 * t14 * t59 - t126)
     #) + t128 * (t2 * t3 * t18 * t54 * t20 + t14 * t31 * t10 * t12 - t7
     #4 * t17 * t54 * t59 + t127 * t59)) - t88 * t14 * t33 * t34 * t22) 
     #+ t5 * (t5 * (t5 * (t5 * (t5 * (t84 * t62 * t33 + t84 * (t2 * (t12
     #1 * t11 + t89) - t122 * t11) * t5) + t121 * t2 * t31 * t12 - t122 
     #* t31 * t12 - t88 * t1 * t15 * t36 + t91 * (t76 + t65)) - t123 * t
     #7 * t21) + t75 * t41 * t33 * (t121 * t54 + t124 * t20)) - t62 * t3
     #4 * t52 * t120) + t93 * t34 * (t14 * (t11 * t120 + t79) + t24 * t1
     #20) + t123 * t54 * t13 * t20
      t79 = t4 * t6
      t84 = t69 * t4
      t91 = t69 * t32
      t93 = t16 * t6 * t11
      t99 = t64 * t1
      t106 = t99 * t3
      t120 = t95 + t94
      t123 = t30 * t5
      t129 = t35 * t53
      t130 = t1 * t31 * t28
      t23 = t111 * (t48 * (t10 * (t3 * (t86 * t3 * t5 * t54 - t55 * t120
     # - t56 * t120) + t110 * t48 * (t123 * t3 - t114)) - t81 * t69 * t6
     #6 * t29) * t9 * t8 - t17 * t11 * (t129 * t7 + t130)) + t106 * (t30
     # * t19 * t4 * t29 * t54 - t115 * t35 * t26 * t22 - t121 * t58 * t2
     #9 + t91 * t17 * t28 + t115 * t36 * (t56 * t11 + t92) * t7 + t128 *
     # t103 * (t79 * t30 * t13 - t16 * t27) + t46 * t2 * t29 * (-t45 + t
     #84)) + t48 * (t5 * (t5 * (t71 * t48 * (t2 * (t23 * t7 + t62) - t93
     #) + t126) + t7 * t6 * t36 * (-t62 * t2 + t24 * t69 + t93)) + t17 *
     # t33 * t53 * t59) * t9 * t8
      t24 = t1 * t36
      t62 = t29 * t7
      t92 = t101 * t48
      t93 = t1 * t32
      t110 = t34 * t53 * t10
      t114 = t31 * t27
      t120 = t37 * t7
      t121 = t10 * t3
      t126 = t115 * t33 * t26
      t32 = t33 * (t125 * t32 * t26 * t8 * t9 * t6 + t82 * t26 * t8 * t9
     # * t59) + t49 * t29 * t30 * (t80 * t59 * t9 * t8 + t65 * t59) * t4
     # + t104 * t64 * t8 * t9 * t59 * t52 + t85 * t10 * t59 * t53 * t36 
     #- t129 * t102 * t67 - t130 * t89 * t10
      t82 = t124 * t26
      t85 = t3 * t33 * t54
      t102 = t121 * t95 * t94
      t104 = t92 * t69 * t41
      t124 = t104 * (t1 * t30 * t29 + t36 * t52 * t7)
      t125 = t104 * (t6 * (t55 * t3 + t128 * (t50 - t80)) - t55 * t66)
      t130 = t55 * t86 * t76 * t11 * t7
      t131 = t49 + t50
      t88 = t64 * t48 * t38 * t6 * (t24 * t88 + t62 * t90)
      t90 = t95 * t98 * t6 * t10 * (t128 * t10 + t108)
      t68 = t109 * t56 * t68 * (t101 * t5 + t74 * t11)
      t72 = t128 * t72 * t64 * (t80 * t65 + t76 * t66)
      t74 = t48 * t37 * t10 * t6 * (t106 * (t85 + t82) + (t127 * t10 + t
     #1 * (t30 * (t48 * t7 + t109) + t56 * t1) * t11 + t101 * (t48 * t1 
     #* t52 + t55 * t7 + t56 * t7) + t113 * t14 * t7) * t9 * t8)
      t83 = t123 * t7 - t83
      t30 = t31 * t18 * t28 * t22 + t129 * t16 * t12 + t100 * (t10 * (t1
     #0 * (-t121 * t31 * t26 - t128 * t42 * t26) + t4 * t3 * (t6 * (t6 *
     # (-t86 * t5 * t7 - t76 * t4) + t103 * t96) - t111 * t30 * t29)) - 
     #t95 * t81 * t26 * t6) * t9 * t8 - t78 * t2 * t33 * t29 * t83 + t77
     # * t36 * t26 * t6 * t83
      t33 = t55 * t56 * (t5 * t18 * t22 + t122 * t12)
      t55 = t126 * (t15 * t2 * t26 * t13 + t128 * t19 * t21)
      t4 = -14336 * t104 * t8 * t9 * t57 - 22528 * t130 * t8 * t9 + 1280
     #0 * t72 + 13824 * t88 + 19968 * t68 + 6656 * t119 + 7424 * t33 + 7
     #68 * t48 * (t6 * (t6 * (t2 * (t3 * (t3 * (t3 * (-t38 * t29 * t20 +
     # t39 * t29 * t20) - t41 * (t16 * t29 + t40)) + t15 * t4 * t29 * t1
     #1) - t43 * t7 * t12) + t44 * t5 * t13 * t6) + t46 * t45 * t5 * t25
     #) + t47 * t19 * t4 * t21) + 1792 * t55 + 2048 * t30 + 5120 * t107 
     #- 8448 * t90 - 8192 * t112 - 7168 * t74 - 11264 * t117 * t8 * t9 *
     # t59 - 9216 * t113 * t73 * t10 * t131 - 8960 * t117 * (t49 * t25 +
     # t50 * t25) - 1024 * t23 - 512 * t51 - 256 * t75 - 1536 * t61 - 40
     #96 * t118 * t9 * t8 * (t115 * (t105 * t25 + t63 * t25) + t102)
      t7 = -2560
      t15 = 3072
      t21 = t80 + t49 - t50 - t66
      t23 = 4
      t25 = 2
      t25 = t25 * (t80 * (t49 + t50 - t66) - t37 * t6 * t10 + t66 * t131
     #)
      t29 = t23 * (t84 * t10 + t120 * t5)
      t30 = t76 + t95 + t65 + t94 + t25 + t29
      t33 = t21 ** 2
      t23 = t23 * (t103 + t79) * (t1 * t10 + t67)
      t21 = -t21 * cdsqrt(t30) + t33 + t23
      t23 = t33 + t23
      t23 = cdsqrt(t23)
      t23 = -t49 * t23 - t80 * t23 - t76 - t66 * (-t23 + t66) + t50 * t2
     #3 - t65 - t94 - t25 - t29
      t25 = 0.1D1 / t9
      t29 = 0.1D1 / t8
      t21 = 0.1D1 / t21 ** 2
      t23 = 0.1D1 / t23 ** 2
      ret = (5376 * t92 * t1 * t11 * (t24 * t53 + t62 * t31) + t15 * (t1
     #9 * (-t91 * t64 * t27 * t8 * t9 + t110 * t87 * t1) + t3 * (-t1 * t
     #35 * t54 * t8 * t9 * t11 * t20 + t71 * t66 * t13 * (-t58 * t8 * t9
     # + t93 * t27)) + t99 * t9 * t8 * (-t114 * t10 - t93 * t28 - t110) 
     #* t17 + t128 * t47 * t16 * t12 + t58 * t2 * t18 * t26 * t22) - 563
     #2 * t62 * t69 * t2 * t36 * (t42 * t13 + t89 * t6) + t7 * (t126 * (
     #t116 * t59 + t60 * t66) * t9 * t8 + t121 * (t6 * (t52 * (t87 * t14
     # * t35 * t10 + t10 * t34 * (-t16 * t5 * t10 + t41 * t14 + t81 * t2
     #0) * t6) + t114 * t3 * t59) + t70 * t27 * (t97 - t120))) - 3584 * 
     #t32 - 3328 * t102 * (t85 + t82) + 7680 * t124 - 16384 * t125 + 176
     #64 * t130 * (t80 + t66) - 17920 * t130 * t131 + t4) * t29 * t25 * 
     #t21 * t23


      hjetmass_qqbqqb_bubble_pmpm_s34_dp =  ret/16d0*(0d0,1d0)
      return
      
      end function

      double complex function
     & hjetmass_qqbQQB_triangle_pmpm_s34_mhsq_s12_dp
     & (i1,i2,i3,i4,za,zb,mt)
      implicit double complex (t)
      integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      double complex ret
      double precision cg
      double precision mt

      cg =  1d0

      t1 = za(i1, i4)
      t2 = za(i2, i3)
      t3 = zb(i3, i1)
      t4 = za(i1, i2)
      t5 = za(i3, i4)
      t6 = zb(i2, i1)
      t7 = zb(i4, i3)
      t8 = za(i1, i3)
      t9 = za(i2, i4)
      t10 = zb(i3, i2)
      t11 = zb(i4, i1)
      t12 = zb(i4, i2)
      t13 = t8 * t3
      t14 = t2 * t10
      t15 = t1 * t11
      t16 = t9 * t12
      t17 = t15 + t16 + t13 + t14
      t17 = -0.4D1 * t4 * t5 * t6 * t7 + t17 ** 2
      t13 = cg * cdsqrt(t17) + t13 + t14 + t15 + t16
      t14 = 0.2D1 * t4 * t6
      t15 = t14 + t13
      t16 = 0.1D1 / 0.2D1
      t14 = t16 * t13 ** 2 - t14 * t5 * t7
      t17 = 0.2D1 * t5 * t7 + t13
      t18 = 0.1D1 / 0.4D1
      t19 = t13 ** 2
      t20 = t4 * t5
      t21 = -t20 * t6 * t7 + t18 * t19
      t21 = 0.1D1 / t21
      t22 = t1 * t3 + t10 * t9
      t19 = t19 * t21
      t23 = t19 * t22
      t24 = t11 * t9 + t2 * t3
      t19 = t19 * t24
      t25 = t13 * t4
      t24 = t25 * t6 * t21 * t24
      t18 = t13 * t18
      t26 = -t16 * t9 * t6 * t7 + t18 * t3
      t27 = t13 * t5 * t21
      t16 = t20 * t16 * t3 - t18 * t9
      t18 = t27 * t7 * t22
      t20 = 0.1D1 / t4
      t22 = 0.1D1 / t5
      t28 = 0.1D1 / t6
      t29 = 0.1D1 / t7
      t14 = 0.1D1 / t14 ** 2
      t30 = t3 ** 2
      t31 = t29 * t28 * t22 * t20
      t32 = t15 * t14 * t17
      t16 = -t13 ** 2 * t6 * t21 ** 2 * t16 ** 2 * t7 - t25 * t21 * t26 
     #** 2 * t27
      t13 = 0.1D1 / t13
      t21 = -0.8D1
      ret = t21 * (t32 * t24 * t18 + t31 * (t9 * (t30 * t8 + t9 * (-t10 
     #* t11 + t12 * t3)) - t1 * t2 * t30)) - 0.32D2 * t32 * t16 - 0.64D2
     # * t31 * mt ** 2 * t13 * t16 - 0.4D1 * t14 * (t5 * t7 * t15 ** 2 *
     # t23 * t24 * t20 * t28 + t4 * t6 * t17 ** 2 * t19 * t18 * t22 * t2
     #9) - 0.2D1 * t32 * t23 * t19

      hjetmass_qqbQQB_triangle_pmpm_s34_mhsq_s12_dp = ret/16d0/(0,1d0)
      return

      end function
     &
      double complex function
     & hjetmass_qqbQQb_triangle_pmpm_s34_mhsq_s12_rat_dp
     & (i1,i2,i3,i4,za,zb)
      implicit double complex (t)
      integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      double complex ret
      double precision cg

      cg = 1d0

      t1 = zb(i3, i1)
      t2 = za(i2, i4)
      t3 = za(i1, i3) * t1
      t4 = za(i2, i3) * zb(i3, i2)
      t5 = za(i1, i4) * zb(i4, i1)
      t6 = t2 * zb(i4, i2)
      t7 = t3 + t4 + t5 + t6
      t8 = za(i1, i2)
      t9 = za(i3, i4)
      t10 = zb(i2, i1)
      t11 = zb(i4, i3)
      t7 = -0.4D1 * t8 * t9 * t10 * t11 + t7 ** 2
      t3 = cg * cdsqrt(t7) + t3 + t4 + t5 + t6
      t4 = t8 * t9
      t5 = t2 ** 2 * t10 * t11 + t4 * t1 ** 2
      t6 = t3 ** 2
      t7 = 0.1D1 / 0.4D1
      t7 = -t4 * t10 * t11 + t6 * t7
      t8 = 0.1D1 / t8
      t12 = 0.1D1 / t11
      t13 = 0.1D1 / t10
      t9 = 0.1D1 / t9
      t7 = 0.1D1 / t7 ** 2
      ret = (-0.32D2 * t4 * t6 * t2 * t10 * t1 * t11 + 0.16D2 * t4 * t3 
     #* t10 * t11 * t5 + 0.4D1 * t3 * t6 * t5) * t8 * t9 * t13 * t12 * t
     #7
    
      hjetmass_qqbQQb_triangle_pmpm_s34_mhsq_s12_rat_dp =
     & ret/16d0/(0,1d0)
      return

      end function

      double complex function
     & hjetmass_qqbQQB_triangle_pmmp_s34_mhsq_s12_dp
     & (i1,i2,i3,i4,za,zb,mt)
      implicit double complex (t)
      integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      double complex ret
      double precision cg
      double precision mt

      cg =  1d0

      t1 = za(i2, i3)
      t2 = zb(i3, i2)
      t3 = zb(i4, i1)
      t4 = za(i1, i2)
      t5 = za(i3, i4)
      t6 = zb(i2, i1)
      t7 = zb(i4, i3)
      t8 = za(i1, i4)
      t9 = za(i1, i3)
      t10 = za(i2, i4)
      t11 = zb(i3, i1)
      t12 = zb(i4, i2)
      t13 = t9 * t11
      t14 = t1 * t2
      t15 = t8 * t3
      t16 = t10 * t12
      t17 = t15 + t16 + t13 + t14
      t17 = -0.4D1 * t4 * t5 * t6 * t7 + t17 ** 2
      t13 = cg * cdsqrt(t17) + t13 + t14 + t15 + t16
      t14 = 0.2D1 * t4 * t6
      t15 = t14 + t13
      t16 = 0.1D1 / 0.2D1
      t14 = t16 * t13 ** 2 - t14 * t5 * t7
      t17 = 0.2D1 * t5 * t7 + t13
      t18 = 0.1D1 / 0.4D1
      t19 = t13 ** 2
      t20 = t4 * t5
      t21 = -t20 * t6 * t7 + t18 * t19
      t21 = 0.1D1 / t21
      t22 = t1 * t11 + t10 * t3
      t19 = t19 * t21
      t23 = t19 * t22
      t24 = t1 * t12 + t3 * t9
      t19 = t19 * t24
      t25 = t13 * t4
      t22 = t25 * t6 * t21 * t22
      t18 = t13 * t18
      t26 = t16 * t1 * t6 * t7 + t18 * t3
      t27 = t13 * t5 * t21
      t16 = t20 * t16 * t3 + t18 * t1
      t18 = t27 * t7 * t24
      t20 = 0.1D1 / t7
      t24 = 0.1D1 / t4
      t28 = 0.1D1 / t6
      t29 = 0.1D1 / t5
      t14 = 0.1D1 / t14 ** 2
      t30 = t1 ** 2
      t31 = t29 * t28 * t24 * t20
      t32 = t15 * t14 * t17
      t16 = t13 ** 2 * t6 * t21 ** 2 * t16 ** 2 * t7 + t25 * t21 * t26 *
     #* 2 * t27
      t13 = 0.1D1 / t13
      t21 = -0.8D1
      ret = t21 * (t32 * t22 * t18 + t31 * (t3 * (t2 * t30 + t3 * (t1 * 
     #t8 - t10 * t9)) - t30 * t11 * t12)) - 0.32D2 * t32 * t16 - 0.64D2 
     #* t31 * mt ** 2 * t13 * t16 - 0.4D1 * t14 * (t5 * t7 * t15 ** 2 * 
     #t22 * t19 * t24 * t28 + t4 * t6 * t17 ** 2 * t23 * t18 * t29 * t20
     #) - 0.2D1 * t32 * t23 * t19

      hjetmass_qqbQQB_triangle_pmmp_s34_mhsq_s12_dp= ret/16d0/(0,1d0)
      return

      end function
     &
      double complex function
     & hjetmass_qqbQQb_triangle_pmmp_s34_mhsq_s12_rat_dp
     & (i1,i2,i3,i4,za,zb)
      implicit double complex (t)
      integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      double complex ret
      double precision cg

      cg = 1d0

      t1 = za(i2, i3)
      t2 = zb(i4, i1)
      t3 = za(i1, i3) * zb(i3, i1)
      t4 = t1 * zb(i3, i2)
      t5 = za(i1, i4) * t2
      t6 = za(i2, i4) * zb(i4, i2)
      t7 = t3 + t4 + t5 + t6
      t8 = za(i1, i2)
      t9 = za(i3, i4)
      t10 = zb(i2, i1)
      t11 = zb(i4, i3)
      t7 = -0.4D1 * t8 * t9 * t10 * t11 + t7 ** 2
      t3 = cg * cdsqrt(t7) + t3 + t4 + t5 + t6
      t4 = t8 * t9
      t5 = t1 ** 2 * t10 * t11 + t4 * t2 ** 2
      t6 = t3 ** 2
      t7 = 0.1D1 / 0.4D1
      t7 = -t4 * t10 * t11 + t6 * t7
      t8 = 0.1D1 / t8
      t12 = 0.1D1 / t11
      t9 = 0.1D1 / t9
      t13 = 0.1D1 / t10
      t7 = 0.1D1 / t7 ** 2
      ret = (-0.32D2 * t4 * t6 * t1 * t10 * t2 * t11 - 0.16D2 * t4 * t3 
     #* t10 * t11 * t5 - 0.4D1 * t3 * t6 * t5) * t8 * t9 * t13 * t12 * t
     #7

      hjetmass_qqbQQb_triangle_pmmp_s34_mhsq_s12_rat_dp =
     & ret/16d0/(0,1d0)
      return

      end function
      double complex function hjetmass_qqbqqb_bubble_pmmp_s34_dp
     &  (i1,i2,i3,i4,za,zb)
      implicit double complex (t)
      integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      double complex ret

      t1 = za(i1, i2)
      t2 = za(i1, i3)
      t3 = za(i2, i3)
      t4 = za(i2, i4)
      t5 = zb(i2, i1)
      t6 = zb(i3, i1)
      t7 = zb(i3, i2)
      t8 = zb(i4, i1)
      t9 = zb(i4, i2)
      t10 = za(i1, i4)
      t11 = t4 ** 2
      t12 = t11 ** 2
      t13 = t4 * t11
      t14 = t8 ** 2
      t15 = t14 ** 2
      t16 = t8 * t15
      t17 = t8 * t14
      t18 = t3 ** 2
      t19 = t18 ** 2
      t20 = t3 * t19
      t21 = t3 * t18
      t22 = t10 ** 2
      t23 = t22 ** 2
      t24 = t10 * t22
      t25 = t7 ** 2
      t26 = t25 ** 2
      t27 = t7 * t25
      t28 = t6 ** 2
      t29 = t28 ** 2
      t30 = t6 * t29
      t31 = t6 * t28
      t32 = t18 * t9
      t33 = t32 * t2
      t34 = t33 * t14
      t35 = t11 * t25
      t36 = t28 * t22
      t37 = t36 + t35
      t38 = t9 ** 2
      t39 = t38 ** 2
      t40 = t9 * t38
      t41 = t2 ** 2
      t42 = t41 ** 2
      t43 = t2 * t41
      t44 = t35 * t2
      t45 = t36 * t2
      t46 = t6 * t8
      t47 = t6 * t19
      t48 = t42 * t4 * t15
      t49 = t2 * t4
      t50 = t49 * t7
      t51 = t8 * t22
      t52 = t43 * t3
      t53 = t52 * t4
      t54 = t10 * t8
      t55 = t4 * t9
      t56 = t55 * t54
      t57 = t2 * t6
      t58 = t3 * t7
      t59 = t58 * t57
      t60 = t6 * t4
      t61 = t11 * t38
      t62 = t41 * t28
      t63 = t41 * t11
      t64 = t43 * t13
      t65 = t18 * t7
      t66 = t35 * t41
      t67 = t65 * t41
      t68 = t28 * t37
      t69 = t18 * t25
      t70 = t69 * t36
      t71 = t18 * t22
      t72 = t3 * t10
      t73 = t10 * t9
      t74 = t73 * t25
      t75 = t61 + t69
      t76 = t13 * t25
      t77 = t3 * t28
      t78 = t21 * t28
      t79 = t78 * t38
      t80 = t22 * t14
      t81 = t43 * t6
      t82 = t81 * t25
      t83 = t4 * t40
      t84 = t83 * t22
      t85 = t36 * t18
      t86 = t38 * t27
      t87 = t47 * t40
      t66 = t8 * (t8 * (t8 * (t8 * (t2 * t42 * t4 * t29 * t22 - t42 * t3
     # * t29 * t24 - t79 * t10 * t23 + t64 * (t62 + t61) * t25 + t80 * t
     #41 * (t2 * (t36 * t4 + t76) - t77 * t24)) - t63 * t19 * t7 * t26) 
     #+ t33 * t11 * t28 * (t84 + t82)) + t86 * t19 * (-t11 * t75 - t85))
     # + t14 * (t8 * (t8 * (t10 * (t21 * (t63 * t26 + t22 * (t62 - t61) 
     #* t25) + t64 * t28 * t38 * t10) + t8 * (-t66 * t3 * t24 * t8 + t67
     # * t22 * t37)) + t65 * (t71 * t37 * t38 - t70 * t41 - t68 * t42)) 
     #+ t32 * (t2 * t4 * t12 * t25 * t40 - t72 * t11 * t40 * t37 + t74 *
     # t37 * t21 + t42 * t30 * t22)) + t66 * t19 * t31 * t40 + t87 * (t4
     #1 * t29 * t22 + t35 * t75 + t36 * t75)
      t75 = t21 * t40
      t88 = t43 * t17
      t89 = t88 - t75
      t90 = t2 * t11
      t91 = t90 * t27
      t92 = t6 * t7
      t93 = t92 * (t57 - t55)
      t94 = t35 * t32
      t95 = t43 * t4
      t96 = t95 * t31
      t97 = t72 * t41
      t98 = t2 * t8
      t99 = t98 * t24
      t100 = t60 * t9
      t101 = t100 * t22
      t102 = t10 * t6
      t103 = t4 * t7
      t104 = t18 * t38
      t105 = t104 * t14
      t106 = t13 * t27
      t107 = t8 * t3
      t108 = t13 * t39
      t109 = t41 * t14
      t110 = t109 * t4
      t111 = t46 * t41
      t112 = t58 * t41
      t113 = t95 * t15
      t114 = t41 * t18
      t115 = t28 * t42
      t116 = t115 * t7
      t117 = t28 * t23
      t118 = t11 * t14
      t119 = t31 * t24
      t114 = t6 * (t10 * (t10 * (t10 * (-t103 * t41 * t14 * t15 * t10 + 
     #t103 * t14 * (t112 * t17 + t81 * t17 + t75 * t6)) + t14 * (t18 * (
     #-t108 * t7 + t110 * t27) - t47 * t2 * t27 * t9 + t111 * t13 * t40)
     #) + t58 * t4 * (t25 * (t83 * t21 * t8 + t113) + t114 * t31 * t40 -
     # t116 * t17 - t114 * t27 * t17)) - t75 * t106 * t98) + t8 * (t41 *
     # (-t118 * t18 * t26 - t117 * t15) + t88 * (t119 - t106) + t75 * (-
     #t119 + t106) + t104 * t50 * t46 * t10 * (t103 - t102)) * t5 * t1
      t120 = t13 * t17
      t121 = t2 * t19
      t117 = t117 * t18 * t17
      t122 = t25 * t17
      t123 = t122 * t12
      t124 = t10 * t14
      t125 = t106 * t9
      t126 = t17 * t22
      t127 = t18 * t28
      t128 = t14 * t3
      t129 = t128 * (t41 * (-t126 * t35 + t107 * (t119 * t9 - t125 + t36
     # * t7 * (-t58 + t55) - t73 * t35 * t6)) + t52 * t46 * t7 * t37 + t
     #127 * t83 * t24) * t5 * t1
      t130 = t6 * t24
      t131 = t46 * t21
      t132 = t2 * t7
      t133 = t25 * t8
      t134 = t21 * t38
      t135 = t90 * t7 * t8
      t70 = t135 * (t41 * (t130 * t15 - t131 * t27) + t124 * t42 * t31 +
     # t32 * (t102 * t40 * t11 - t21 * t26)) + t134 * (-t99 * t31 - t80 
     #* t35 - t70 + t133 * (t132 + t73) * t13) * t5 * t1
      t30 = t8 * t70 + t9 * (t9 * (t9 * (t3 * (t7 * (t7 * (t103 * t19 * 
     #t28 * (t72 + t49) + t120 * (-t71 + t63)) + t121 * t29 * t22) - t11
     #7 * t4) + t19 * t13 * t28 * t7 * t10 * t38) + t51 * t43 * t21 * t3
     #0) + t124 * t2 * (t41 * (t36 * t17 * t11 + t123) - t117 + t95 * t6
     #5 * t29)) + t3 * t114 + t8 * (t8 * (t8 * (t8 * (t99 * (-t97 * t31 
     #- t94 + t96) + t73 * (t21 * (t91 * t10 + t93 * t24) + t81 * t12 * 
     #t7 * t9 - t53 * t31 * t22)) + t74 * t19 * (t44 + t45 + t101)) + t8
     #6 * t60 * t20 * t22) + t38 * t25 * t4 * t21 * (-t102 * t21 * t25 -
     # t94 + t96)) + t106 * t20 * t6 * t39 + t107 * (t28 * (t103 * t89 *
     # t22 - t105 * t23) - t61 * t19 * t26 - t102 * t35 * t89) * t5 * t1
     # - t88 * t22 * t21 * t31 * t25 + t129
      t70 = t65 * t38
      t71 = t81 * t14
      t81 = t83 * t18
      t83 = t60 * t10 * t7
      t86 = t3 * t9
      t89 = t90 * t1 * t5
      t68 = t18 * (t68 * t1 * t43 * t5 * t14 * t9 + t89 * t14 * t37 * t4
     #0) - t76 * t102 * t42 * t16 - t103 * t22 * t20 * t31 * t39 + t55 *
     # t17 * t41 * (t54 * t37 * t5 * t1 + t68 * t41) * t3 + t131 * t90 *
     # t37 * t39 + t121 * t46 * t1 * t5 * t7 * t37 * t38
      t99 = t17 * t4
      t114 = t36 + t35
      t117 = t5 * t1
      t101 = t122 * t101 * t2
      t121 = t117 * t7
      t33 = t54 * (t19 * (t73 * t60 * t2 * t27 * t14 + t121 * t8 * t114 
     #* t38) + t21 * (t117 * t11 * t28 * t10 * t39 + t95 * t29 * t7 * t3
     #8 + t101) - t116 * t10 * t11 * t15 + t112 * t12 * t6 * t14 * t40 -
     # t130 * t33 * t17 * (t103 * t8 + t117 * t6)) + t83 * (-t2 * t20 * 
     #t27 * t14 * t9 + t73 * t43 * t11 * t16 - t82 * t21 * t17 - t84 * t
     #21 * t17 + t93 * t20 * t40) + t107 * (t7 * (t7 * (t2 * (-t36 * t21
     # * t9 * t8 - t86 * t126 * t11) + t7 * (-t135 * t21 * t9 + t72 * t4
     # * (t41 * (-t128 * t6 + t99) - t134 * t6)) + t18 * t12 * t39 + t11
     #8 * t115 + t124 * t53 * t28) + t46 * t22 * (-t99 * t41 * t10 - t10
     #4 * t54 * t4 + t18 * t11 * t40 + t97 * t6 * t14)) + t80 * t42 * t2
     #9) * t5 * t1
      t82 = -t41 * t7 * t8 + t72 * t38
      t84 = t73 - t132
      t93 = t125 * t21 * t14
      t82 = -t119 * t2 * t18 * t17 * t82 - t42 * t12 * t27 * t16 - t24 *
     # t20 * t29 * t39 + t93 * t82 + t117 * t58 * t54 * (t8 * (t8 * (t2 
     #* (t32 * t37 + t96) - t72 * t49 * t46 * t84) + t100 * t21 * t7 * t
     #84) + t108 * t18 * t6)
      t84 = t7 * t40
      t43 = t31 * (t117 * t50 * t19 * t40 * t10 - t58 * t48 * t22) + t6 
     #* (t121 * t43 * t11 * t16 * t22 + t54 * t13 * t7 * t9 * (t117 * t8
     #8 - t84 * t19)) - t75 * t2 * t12 * t27 * t14 - t88 * t32 * t24 * t
     #29 + t117 * t103 * (t19 * t4 * t39 + t84 * t20 + t48) * t10 * t28
      t75 = t127 + t118
      t84 = t80 * t69 * t60
      t88 = t90 * t28
      t52 = t107 * (t9 * (t9 * (t2 * (t10 * (t10 * (t10 * (t99 * t72 * t
     #28 - t88 * t17) + t58 * (-t2 * t18 * t29 + t120 * t7)) - t123 * t2
     #) - t91 * t78) + t21 * t31 * t8 * t23 * t9) + t111 * t21 * t25 * t
     #114) + t120 * t52 * t26)
      t19 = t19 * t31 * t40
      t78 = t7 * t11
      t90 = t117 * t124
      t91 = t55 + t57
      t96 = t35 * t91 + t36 * t91
      t97 = t58 * t54
      t16 = t117 * (t95 * t10 * t37 * t16 + t92 * t40 * t37 * t20 + t113
     # * t96 + t87 * t96) + t97 * (t2 * (t61 * t46 * t3 * (t69 - t80) + 
     #t85 * t14 * t38 * t4) + t41 * (t120 * t102 * t38 + t94 * t6 * t14 
     #+ t77 * t55 * (t18 * (t92 * t9 - t133) + t126)) + t19 * t22 + t64 
     #* t25 * t15)
      t18 = t55 + t57
      t40 = t32 * t6
      t64 = t84 * (t81 + t71)
      t39 = t97 * t11 * t28 * (t21 * t10 * t39 + t42 * t7 * t17)
      t42 = t97 * t89 * t28 * t9 * (t109 + t104)
      t77 = t72 * t50 * t46 * t9 * (t63 * t7 * t17 + t79 * t10)
      t79 = t97 * t88 * t9
      t85 = t79 * (t9 * (t32 * (-t54 + t57) + t110) - t112 * t14)
      t79 = t79 * (t41 * t10 * t17 + t134 * t7)
      t87 = t70 * t124 * t63 * t28
      t88 = t101 * t21 * t91
      t89 = t105 * t41 * (t8 * t12 * t27 + t3 * t24 * t29)
      t67 = t124 * t67 * t60 * t38 * (t72 * t28 + t78 * t8)
      t92 = t60 * t34 * t22 * t25 * (t49 * t14 + t40)
      t94 = t58 + t54
      t4 = 54 * t77 + 69 * t87 * t94 + 78 * t67 + 29 * t89 + 30 * t79 + 
     #50 * t83 * t34 * (t62 * t54 + t61 * t58) + 21 * t39 + 26 * t86 * t
     #98 * (t119 * t104 * t4 * t8 + t109 * t106 * t3 * t6 + t41 * t12 * 
     #t27 * t17 + t134 * t24 * t29) + 3 * t52 - 6 * t16 - 13 * t64 - 22 
     #* t132 * t73 * t21 * t17 * (t2 * t22 * t31 + t76 * t9) - 16 * t117
     # * t86 * t98 * (t84 + (t35 * t75 + t36 * t75) * t9 * t2) - 33 * t9
     #2 - 32 * t111 * t38 * t4 * t3 * (t121 * t10 * t75 + t49 * t114 * t
     #14 + t40 * t114) - 56 * t42 - 44 * t117 * t105 * t60 * t41 * t37 -
     # 36 * t88 - 35 * t105 * t60 * t41 * (t35 * t18 + t36 * t18) - 70 *
     # t87 * t91 - 64 * t85 - 88 * t90 * t70 * t63 * t28
      t12 = -10
      t16 = t57 + t58 - t54 - t55
      t17 = 4 * t72 * t6 * t9 + 4 * t50 * t8
      t18 = 2 * t55 * t94 - 2 * t97 + 2 * t57 * (t58 + t54 - t55)
      t27 = t17 + t18 + t62 + t61 + t69 + t80
      t29 = t16 ** 2
      t35 = (4 * t102 + 4 * t103) * (t86 + t98)
      t16 = t16 * cdsqrt(t27) + t29 + t35
      t27 = t29 + t35
      t27 = cdsqrt(t27)
      t17 = -t17 - t18 - t62 - t61 + t58 * t27 + t57 * t27 - t54 * (t54 
     #+ t27) - t55 * t27 - t69
      t18 = 0.1D1 / t5
      t27 = 0.1D1 / t1
      t29 = 256
      t16 = 0.1D1 / t16 ** 2
      t17 = 0.1D1 / t17 ** 2
      ret = t29 * (7 * t34 * (t2 * t23 * t31 * t14 + t32 * t13 * t26) + 
     #t12 * (t60 * (t9 * (t38 * (t46 * t20 * t25 * t22 + t47 * (-t46 * t
     #24 + t44 + t45) * t9) + t48 * t37) + t53 * t25 * t15 * (t51 - t50)
     #) + t34 * (t56 * t37 + t59 * t37) * t5 * t1) - 2 * t30 - 28 * t86 
     #* t49 * t46 * (t83 * (t81 + t71) + (t74 * t21 * t6 + t10 * (t41 * 
     #(t58 * t8 + t124) + t104 * t10) * t28 + t7 * (t54 * t3 * t38 + t10
     #9 * t7 + t70) * t11 + t126 * t50) * t5 * t1) - 14 * t68 - 4 * t33 
     #- 8 * t82 - 12 * t43 - 20 * t9 * t2 * (t10 * (t78 * (t69 * t49 * t
     #15 + t113 * t28 + t19) + t80 * t32 * t31 * (-t49 * t8 + t65)) - t9
     #3 * t57 + t90 * t65 * t60 * (t59 + t56)) - t66 + t4) * t27 * t18 *
     # t16 * t17

      hjetmass_qqbqqb_bubble_pmmp_s34_dp = ret/16d0*(0d0,1d0)
      return
      
      end function

      double complex function hjetmass_qqbqqb_bubble_pmmp_s12_dp
     &  (i1,i2,i3,i4,za,zb)
      implicit double complex (t)
      integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      double complex ret

      t1 = za(i1, i3)
      t2 = za(i1, i4)
      t3 = za(i2, i3)
      t4 = za(i2, i4)
      t5 = zb(i3, i1)
      t6 = zb(i4, i1)
      t7 = zb(i4, i2)
      t8 = za(i3, i4)
      t9 = zb(i3, i2)
      t10 = zb(i4, i3)
      t11 = t5 ** 2
      t12 = t11 ** 2
      t13 = t5 * t11
      t14 = t1 ** 2
      t15 = t14 ** 2
      t16 = t1 * t15
      t17 = t1 * t14
      t18 = t3 ** 2
      t19 = t18 ** 2
      t20 = t3 * t19
      t21 = t3 * t18
      t22 = t9 ** 2
      t23 = t22 ** 2
      t24 = t9 * t22
      t25 = t14 * t11
      t26 = t18 * t22
      t27 = t25 + t26
      t28 = t6 ** 2
      t29 = t28 ** 2
      t30 = t6 * t29
      t31 = t6 * t28
      t32 = t2 ** 2
      t33 = t32 ** 2
      t34 = t2 * t32
      t35 = t4 ** 2
      t36 = t35 ** 2
      t37 = t4 * t35
      t38 = t7 ** 2
      t39 = t38 ** 2
      t40 = t7 * t38
      t41 = t1 * t4
      t42 = t2 * t3
      t43 = t38 * t32
      t44 = t43 * t36
      t45 = t17 * t4
      t46 = t45 * t22
      t47 = t34 * t3
      t48 = t34 * t21
      t49 = t17 * t37
      t50 = t35 * t28
      t51 = t37 * t7
      t52 = t51 * t22
      t53 = t1 * t13
      t54 = t53 * t32
      t55 = t14 * t5
      t56 = t18 * t38
      t57 = t14 * t22
      t58 = t7 * t13
      t59 = t58 * t19
      t60 = t57 + t43
      t61 = t43 * t22
      t62 = t32 * t18
      t63 = t35 * t2
      t64 = t38 * t60
      t65 = t18 * t11
      t66 = t18 * t9
      t67 = t2 * t5
      t68 = t67 * t22
      t69 = t42 * t14
      t70 = t5 * t18
      t27 = t6 * (t6 * (t6 * (t6 * (t6 * (-t47 * t57 * t35 * t6 + t62 * 
     #t35 * t9 * t60) + t63 * (t21 * (t14 * t23 + t61) + t45 * t2 * t11 
     #* t38)) + t66 * (-t61 * t18 * t35 + t65 * t32 * t60 - t64 * t36)) 
     #+ t70 * (t16 * t4 * t13 * t22 + t32 * t36 * t7 * t39 - t69 * t13 *
     # t60 + t68 * t60 * t21)) - t43 * t18 * t19 * t11 * t24) + t6 * (t6
     # * (t6 * (t6 * (t38 * (t17 * t4 * t36 * t22 + t44 * (-t42 + t41) -
     # t2 * t33 * t21 * t11) + t25 * t22 * (-t48 + t49) + t50 * t32 * (-
     #t47 * t38 + t41 * t43 + t46)) - t14 * t19 * t35 * t9 * t23) + t56 
     #* t55 * t4 * (t54 + t52)) - t25 * t19 * t24 * t27) + t59 * (t32 * 
     #(t27 * t38 + t35 * t39) + t57 * t27) + t57 * t19 * t35 * t13 * t40
      t61 = t57 + t43
      t71 = t14 * t28
      t72 = t70 * t7
      t73 = t2 * t9
      t74 = t3 * t11
      t75 = t74 * t7 * t6
      t76 = t17 * t5
      t77 = t32 * t4
      t78 = t4 * t9
      t79 = t1 * t5
      t80 = t4 * t7
      t81 = t80 - t79
      t82 = t2 * t6
      t83 = t42 * t11
      t84 = t14 * t4
      t85 = t32 * t31
      t86 = t37 * t31
      t87 = t86 * t1
      t88 = t65 * t7
      t89 = t50 * t7
      t90 = t41 * t5
      t91 = t90 * t28
      t92 = t6 * t18
      t93 = t1 * t81
      t94 = t56 * t5 * t28
      t95 = t25 * t3
      t96 = t4 * t6
      t97 = t15 * t37
      t98 = t42 * t7
      t99 = t53 * t2
      t100 = t53 * t19
      t101 = t2 * t7
      t102 = t17 * t13
      t103 = t65 + t50
      t104 = t3 * t5
      t105 = t104 * t96
      t106 = t2 * t28
      t107 = t35 * t31
      t108 = t21 * t13
      t109 = t3 * t9
      t110 = t34 * t38
      t111 = t41 * t28
      t112 = t6 * t3
      t113 = t1 * t9
      t114 = t51 * t31
      t115 = t95 * t9
      t116 = t34 * t31
      t117 = t19 * t1
      t118 = t112 * (t14 * (-t18 * t103 * t23 + t101 * t96 * (-t50 + t65
     # - t105) * t22) + t109 * t79 * t43 * (t50 - t65 - t105) + (t108 - 
     #t86) * t24 * t17 + t110 * (-t106 * t65 - t58 * t21 + t107 * (t80 -
     # t82))) * t10 * t8
      t119 = t78 * t55
      t120 = t37 * t40
      t121 = t120 * t32
      t122 = t32 * t28
      t123 = t35 * t38
      t124 = t14 * t37 * t29
      t16 = t9 * (t31 * (t16 * t3 * t35 * t13 - t102 * t32 * t21 + t97 *
     # t67 * t28 - t121 * t21) + t9 * (t9 * (-t114 * t14 * t21 - t119 * 
     #t20 * t28) + t53 * t20 * t38 * (t42 + t41))) + t101 * (t17 * t19 *
     # t5 * t12 * t7 + t2 * t20 * t4 * t13 * t40 + t124 * (t122 + t123))
      t125 = t34 * t18 * t4
      t126 = t98 * t6 * (t32 * (t92 * t58 * t1 - t96 * t56 * t11) + t73 
     #* (-t19 * t11 * t9 - t66 * t50 + t87) * t7 + t108 * t57) * t10 * t
     #8
      t127 = t63 * t31
      t128 = t96 * t5
      t129 = t35 * t6
      t130 = t129 * t7
      t131 = t78 * t3
      t132 = t109 * t67 * t41 * t7
      t74 = t28 * t3 * ((t9 * (t9 * (t14 * t2 * (t18 * (-t82 * t11 + t53
     #) - t127) + t131 * t14 * (t1 * (-t128 + t74) + t130)) + t121 * t11
     #2) + t129 * t47 * t5 * t40) * t10 * t8 + t132 * (t120 + t102))
      t16 = t74 + t9 * t16 + t7 * (t7 * (t113 * t47 * t37 * t30 + t114 *
     # t2 * (t32 * (-t112 * t79 + t111) - t47 * t28 - t26 * t41)) + t117
     # * t11 * t22 * (t116 + t115)) + t9 * (t9 * (t9 * (t73 * t6 * t19 *
     # t1 * (-t89 + t91 - t88) + t92 * (t7 * (t1 * (t3 * (t3 * (-t84 * t
     #13 + t83 * (t79 + t82)) + t85 * t35) + t87 * t2) - t77 * t21 * t5 
     #* t6 * t7) + t77 * t55 * t3 * t31)) + t96 * t18 * (t34 * (t93 * t2
     #9 + t94) + t95 * t35 * t40)) + t101 * (t28 * (t28 * (t5 * (-t48 * 
     #t79 + t80 * t48 + t97 * t5) - t47 * t50 * t1) + t99 * t21 * (t98 -
     # t55)) + t100 * t35 * t40)) + t102 * t43 * t3 * t35 * t31 + t118 -
     # t102 * t20 * t23 * t6 + t126 + t43 * t6 * t5 * (-t122 * t1 * t21 
     #* t11 + t120 * t21 * t5 + t124 * t2 - t125 * t29)
      t55 = t42 * t1
      t74 = t15 * t4
      t87 = t21 * t24
      t95 = t2 * t38
      t97 = t131 * t7
      t98 = t67 * t1
      t102 = t98 * t6
      t114 = t5 * t4
      t118 = t51 * t28
      t53 = t53 * t18
      t120 = t62 * t1 * t22 * t28 * t7
      t121 = t88 * t50 * t1
      t124 = t1 * t36
      t126 = t124 * t29
      t131 = t13 * t19
      t51 = t2 * (t76 * t6 * t9 * (-t86 * t8 * t10 + t131 * t9) * t7 - t
     #100 * t78 * t8 * t10 * t40 + t113 * t10 * t8 * (-t20 * t13 * t9 - 
     #t117 * t12 - t126) * t38) + t32 * (-t51 * t14 * t30 * t8 * t9 * t1
     #0 + t126 * t109 * t40) + t74 * t108 * t24 * t28 + t86 * t70 * t34 
     #* t39
      t100 = t7 * t11
      t64 = t36 * (-t101 * t17 * t22 * t30 + t79 * t64 * t3 * t31) + t4 
     #* (t71 * t18 * t13 * t8 * t10 * t60 + t100 * t21 * (t109 * t60 * t
     #10 * t8 + t25 * t60) * t6) + t98 * t3 * t8 * t10 * t60 * t29 * t35
     # - t113 * t32 * t20 * t12 * t40 + t94 * t8 * t10 * t60 * t37
      t94 = t3 * t1
      t98 = t67 * t7
      t108 = t70 * t4
      t117 = t108 * t28
      t126 = t7 * t1
      t133 = t126 * (t6 * (t31 * (t112 * t1 * t32 * t37 * t22 + t124 * (
     #-t94 * t24 + t43 * t5 + t57 * t5)) + t59 * t32 * (t22 * t3 - t98))
     # + t80 * t19 * t12 * t60) + t117 * (t102 * t60 + t97 * t60) * t10 
     #* t8
      t48 = t48 * t11
      t134 = t122 + t26
      t135 = t14 * t12
      t136 = t58 * t20
      t126 = t126 * t73
      t137 = t109 * t6
      t13 = t137 * (t112 * (t14 * (-t108 * t24 + t2 * (t65 + t50) * t22)
     # + t110 * t50) * t10 * t8 + t98 * t41 * (t45 * t11 * t28 + t70 * t
     #35 * t40 + t112 * (t92 * t2 * t22 + t122 * t109 - t116 - t87))) + 
     #t126 * (-t93 * t37 * t30 * t2 + t136 * t81 * t9 - t52 * t21 * t31 
     #- t54 * t21 * t31) + t112 * (t7 * (t7 * (t32 * (-t128 * t22 * t21 
     #+ t135 * t18) + t55 * t37 * t28 * t22 + t57 * t36 * t28 + t44 * t2
     #8 + t48 * t6 * t9 - t104 * t33 * t4 * t31) + t73 * t1 * (-t50 * t1
     #34 - t65 * t134 + t92 * t99)) + t104 * t57 * (t14 * t3 * t13 - t77
     # * t31)) * t10 * t8
      t44 = t69 * t9 * t6 * t38 * (t2 * t21 * t12 + t31 * t9 * t36)
      t45 = t79 * t60 + t80 * t60
      t52 = t1 * t37
      t54 = t42 * t9 * t6
      t32 = (t52 * t2 * t60 * t30 + t136 * t9 * t60 + t52 * t45 * t29 + 
     #t59 * t45) * t10 * t8 + t54 * (t1 * (t65 * t43 * t28 * t4 + t123 *
     # t104 * (t18 * (t7 * t5 * t9 - t22 * t6) + t85)) + t14 * (t75 * (-
     #t122 + t26) * t4 + t72 * t50 * t22) + t17 * (t29 * t22 * t37 + t10
     #0 * t127) + t131 * t32 * t40)
      t37 = t65 * t2
      t43 = t129 * t9 - t83
      t45 = -t78 + t67
      t12 = t34 * t20 * t12 * t39 + t15 * t36 * t24 * t30 - t125 * t31 *
     # t43 * t40 + t87 * t76 * t28 * t43 + t54 * (t7 * (t7 * (-t114 * t6
     #2 * t28 - t124 * t28 * t7) + t94 * (t6 * (t106 * t4 * t45 - t66 * 
     #t5 * t45) - t135 * t3)) - t117 * t57) * t10 * t8
      t20 = t105 * t10 * t8 * (t114 * (t56 * t61 + t71 * t61) + t120)
      t30 = t119 * t42
      t36 = t65 * t50 * (t15 * t24 * t6 + t47 * t39)
      t43 = t117 * (t4 * t28 * t33 * t40 + t70 * t17 * t23)
      t45 = t14 * t9
      t47 = t80 + t79
      t52 = t73 * t56 * t50 * t25
      t57 = t90 * t85 * t21 * t22 * t7 * t47
      t59 = t91 * t62 * t22 * t7 * (t72 + t111)
      t62 = t132 * t6 * (t95 * t21 * t11 + t45 * t107)
      t63 = t98 * t111 * t66 * (t63 * t6 * t38 + t115)
      t69 = t30 * t6 * t38
      t81 = t69 * (t6 * (t129 * (t109 - t79) + t37) - t88 * t4)
      t88 = t109 + t82
      t11 = -70 * t52 * t47 - 56 * t30 * t8 * t6 * t38 * t10 * t103 - 44
     # * t121 * t8 * t10 * t60 - 36 * t57 - 33 * t59 - 32 * t75 * t35 * 
     #t1 * (t41 * t61 * t28 + t72 * t61 + t73 * (t71 + t56) * t10 * t8) 
     #- 28 * t112 * t90 * t7 * ((t68 * t21 * t7 + t85 * t41 * t9 + t95 *
     # (t35 * (t137 + t106) + t37) + t9 * (t50 * t9 + t83 * t6 + t65 * t
     #9) * t14) * t10 * t8 + t126 * (t53 + t118)) - 22 * t78 * t67 * t21
     # * t31 * (t76 * t22 + t77 * t40) - 16 * t20 - 6 * t32 - 2 * t16 + 
     #7 * t43 + 8 * t12 + 29 * t36 + 30 * t69 * (t21 * t11 * t9 + t127) 
     #+ 50 * t63 + 54 * t62 + 64 * t81 + 69 * t52 * t88 + 78 * t121 * t7
     #3 * (t42 * t38 + t45 * t6) - 88 * t52 * t8 * t10 - t27
      t12 = -t80 + t79 - t109 + t82
      t16 = 4 * t41 * t9 * t6 + 4 * t98 * t3
      t20 = 2 * t80 * t88 - 2 * t54 + 2 * t79 * (-t80 + t109 + t82)
      t27 = t25 + t122 + t123 + t26 + t16 + t20
      t30 = t12 ** 2
      t32 = (4 * t104 + 4 * t96) * (t101 + t113)
      t12 = t12 * cdsqrt(t27) + t30 + t32
      t27 = t30 + t32
      t27 = cdsqrt(t27)
      t16 = -t16 - t20 - t109 * t27 + t79 * t27 + t82 * t27 - t122 - t25
     # - t26 - t80 * (t27 + t80)
      t20 = 0.1D1 / t10
      t25 = 0.1D1 / t8
      t16 = 0.1D1 / t16 ** 2
      t12 = 0.1D1 / t12 ** 2
      ret = (-768 * t112 * (t5 * (t5 * (t4 * (t2 * (t2 * (t2 * (-t55 * t
     #31 * t38 + t84 * t31 * t38) + t109 * (-t17 * t31 * t9 + t18 * t4 *
     # t39)) + t74 * t31 * t22) + t87 * t84 * t38) - t33 * t21 * t6 * t4
     #0 * t5) - t130 * t21 * t22 * t61) - t86 * t17 * t3 * t23) - 5120 *
     # t114 * (t9 * (t95 * (t19 * (t122 * t7 * t5 + t58 * t14) + t49 * t
     #29) + t46 * t18 * t28 * (-t104 * t7 + t106)) - t116 * t70 * t41 * 
     #t40 + t101 * t66 * t1 * t28 * (t102 + t97) * t10 * t8) - 3328 * t1
     #20 * (t53 + t118) - 8960 * t121 * (t79 * t61 + t80 * t61) + 3072 *
     # t51 - 3584 * t64 - 2560 * t133 + 6656 * t105 * (t65 * t1 * t34 * 
     #t6 * t40 + t89 * t17 * t3 * t24 + t107 * t15 * t24 + t48 * t39) - 
     #1024 * t13 + 5376 * t44 + 256 * t11) * t25 * t12 * t16 * t20

      hjetmass_qqbqqb_bubble_pmmp_s12_dp = ret/16d0*(0d0,1d0)
      return
      
      end function

      double complex function hjetmass_qqbqqb_bubble_pmmp_mhsq_dp
     &  (i1,i2,i3,i4,za,zb)
      implicit double complex (t)
      integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      double complex ret
      double complex hjetmass_qqbqqb_bubble_pmmp_s12_dp
      double complex hjetmass_qqbqqb_bubble_pmmp_s34_dp

      hjetmass_qqbqqb_bubble_pmmp_mhsq_dp =
     &   - hjetmass_qqbqqb_bubble_pmmp_s12_dp(i1,i2,i3,i4,za,zb)
     &   - hjetmass_qqbqqb_bubble_pmmp_s34_dp(i1,i2,i3,i4,za,zb)
      return
      
      end function


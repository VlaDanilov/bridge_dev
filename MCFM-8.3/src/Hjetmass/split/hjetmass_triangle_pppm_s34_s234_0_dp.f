
      double complex function hjetmass_triangle_pppm_s34_s234_0_dp 
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision mt
          double precision cg

      t1 = za(i3, i4)
      t2 = za(i2, i3)
      t3 = zb(i3, i2)
      t4 = za(i2, i4)
      t5 = zb(i4, i2)
      t6 = t2 * t3
      t7 = t4 * t5
      t8 = t7 + t6
      if ( dreal(t8) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t8 = cg * cdsqrt(t8 ** 2) + t6 + t7
      t9 = zb(i2, i1)
      t10 = za(i1, i4)
      t11 = zb(i4, i3)
      t12 = za(i1, i2)
      t13 = 0.1D1 / t8
      t14 = t1 * t12
      t15 = 2 * t14
      t16 = t11 * (t15 * t3 * t13 + t10)
      t17 = za(i1, i3)
      t18 = zb(i3, i1)
      t19 = zb(i4, i1)
      t20 = t17 * t18
      t21 = t10 * t19
      t22 = -t15 * t9 * t11 * t13 + t20 + t21
      t23 = t10 * t5 + t17 * t3
      t24 = t12 * t22
      t15 = (-t16 * t2 - t4 * t11 * (t15 * t5 * t13 - t17) + t24) * t9
      t25 = t1 * (2 * t4 * t9 * t11 * t13 - t18)
      t26 = 2 * t2
      t27 = t1 * (t26 * t9 * t11 * t13 + t19)
      t28 = t9 * t22
      t29 = (t25 * t5 + t27 * t3 - t28) * t12
      t30 = t1 * t11
      t31 = 0.1D1 / t12
      t8 = 0.1D1 / t8
      t32 = 0.1D1 / t10
      t33 = 0.1D1 / t22
      t15 = 0.1D1 / t15
      t2 = 0.1D1 / t2
      t34 = t9 * t12
      t35 = (t34 + t7) * t2
      t36 = t35 + t3
      t37 = t15 * t4
      t36 = t37 * t36 * t32 * t23 + t15 * t36 * t22 - t2
      t35 = -t35 - t3
      t38 = t33 ** 2
      t39 = t9 ** 2
      t40 = t9 * t39
      t41 = t4 ** 2
      t42 = t4 * t41
      t43 = t3 ** 2
      t44 = t9 * t16
      t45 = t23 * t18
      t46 = t3 * t22
      t47 = t39 * t16 ** 2 * t38
      t48 = t1 * t32
      t49 = t48 * t2
      t50 = t15 * t39
      t51 = 0.1D1 / t11
      t52 = t32 * (t20 + t30) + t19
      t53 = mt ** 2
      t54 = t8 ** 2
      t55 = t30 * t54
      t56 = t53 * t13 ** 2
      t57 = t56 + t55
      t58 = t15 ** 2
      t59 = t15 * t58
      t24 = t24 * t2
      t60 = t3 * t32
      t61 = t60 * t2
      t62 = t46 * t32
      t28 = t28 * t15
      t63 = t55 * t42 * t40
      t64 = t2 * t58
      t65 = t64 * t39
      t6 = t41 * (t65 * t5 * (-t62 * t57 + t44 * t8 * (t52 * t15 * t22 +
     # t30 * t8 * t32)) * t23 * t41 + t9 * (t39 * (t60 * t55 * (-t24 + t
     #16) * t58 + t46 * t8 * t16 * (t32 * (t30 * (-t6 * t8 + 1) + t20) +
     # t19) * t59) + t40 * (t14 * t11 * t16 * t54 * t32 * t2 * t58 + t24
     # * t8 * t16 * t52 * t59) + t9 * (-t55 * t43 * t22 * t32 * t58 + t6
     #1 * t57 * t15) - t55 * t12 ** 2 * t39 ** 2 * t22 * t16 * t32 * t2 
     #* t59 - t56 * t61 * t31 * t33) * t23 * t4 - t63 * t5 ** 2 * t23 * 
     #t22 * t16 * t32 * t2 * t59 + t53 * t43 * t13 * t32 * t2 * t51 * (-
     #t28 + t31))
      t11 = 0.1D1 / t29
      t1 = 0.1D1 / t1
      t24 = 0.1D1 / t27
      t25 = t25 ** 2
      t27 = t48 * t3
      t29 = t39 * t16 * t22 * t58
      t40 = t9 * t1
      t52 = t51 * t32
      t5 = t52 * (t42 * (t40 * t23 * (t16 * t31 ** 2 * t33 - t29 + t46 *
     # t9 * t30 * (-t26 * t3 * t13 + 1) * t58) * t2 + t3 * t5 * (-t28 + 
     #t31) * t2 ** 2) - t43 * t25 * t24 * (t34 * t11 + t33)) * t13 * t53
     # + t15 * t23 * t9 * t41 * (t39 * (t2 * t8 * t16 * (-t19 * t4 + t32
     # * (t18 * (-t17 * t4 - t14) - t30 * t4)) * t15 + t64 * t22 * t16 *
     # t4 * (t20 * (-t19 * t51 * t1 - t32) - t19)) + t8 * (-t27 * t16 * 
     #t18 + (t46 * (t32 * (t20 + t30) + t19) - t48 * t5 * t16 * t18) * t
     #2 * t4) * t15 * t9 - t49 * t8 * t18 * t3)
      t13 = t46 + t44 - t45
      t14 = t46 - t45
      t26 = t16 * t31
      t28 = t26 * t33
      t53 = t3 * t33
      t54 = t18 * t51
      t55 = t44 - t45
      t11 = t9 * (t12 * t3 * t25 * t51 * t24 * t11 * (t40 - t60) + t54 *
     # (-t46 * t15 + t28) * t2 * t4 + t51 * (t32 * (-t43 * t22 * t15 + t
     #26 * (-t45 * t38 + t53)) + t40 * (t15 * (t20 * t22 * t14 * t15 + t
     #21 * t22 * t15 * t14 - t45 + t46) + t28) + t29 * (t21 + t20) * t1)
     # * t2 * t41 + t64 * t23 * t9 * ((t20 * t13 * t32 + t13 * t19) * t1
     # * t51 + t62) * t42) + t53 * t25 * t51 * t24 * (t40 - t60) + t54 *
     # (-t50 * t16 + t3 * t31) * t2 * t4 + (t3 * t39 * t22 ** 2 * t58 + 
     #t39 * t55 * t58 * t22 + t52 * t31 * (t47 + t43)) * t2 * t41 + t65 
     #* t23 * t32 * t55 * t42
      ret = -8 * t8 * t41 * (t50 * (t44 * t36 + t45 * (t37 * t35 * t32 *
     # t23 + t35 * t15 * t22 + t2) + t46 * t36) + t49 * (t45 * t33 * (t4
     #4 * t33 + t3) - t43 - t47) * t31) - 64 * t6 - 32 * t5 + 4 * t11 + 
     #128 * t63 * t59 * t32 * t23 * t22 * t16 * (t34 * t3 + t7 * (t34 * 
     #t2 + t3)) + 16 * t44 * t2 * t41 * (t27 * t8 * t31 * t33 + t45 * (t
     #51 * (-t20 * t32 - t19) - t48) * t58 * t9 + t22 * t4 * (t30 * t32 
     #+ (t17 ** 2 * t18 ** 2 * t32 + t10 * t19 ** 2) * t1 * t51) * t59 *
     # t23 * t39)

      hjetmass_triangle_pppm_s34_s234_0_dp = ret/32d0/(0,1d0)
      return

      end function

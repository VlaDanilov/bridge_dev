
      double complex function hjetmass_triangle_ppmm_s23_s123_0_dp 
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

      t1 = za(i1, i3)
      t2 = za(i3, i4)
      t3 = zb(i2, i1)
      t4 = zb(i3, i1)
      t5 = za(i1, i4)
      t6 = zb(i4, i1)
      t7 = za(i1, i2)
      t8 = za(i2, i3)
      t9 = zb(i3, i2)
      t10 = t7 * t3
      t11 = t1 * t4
      t12 = t10 + t11
      if ( dreal(t12) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t12 = cg * cdsqrt(t12 ** 2) + t10 + t11
      t13 = 0.1D1 / t12
      t14 = -2 * t10 * t13
      t15 = t8 * t9
      t16 = t15 * (1 + t14)
      t17 = zb(i4, i3)
      t18 = t8 * (-2 * t13 * t6 * t7 * t9 - t17)
      t14 = t11 * t14 * t15 + t3 * (t16 * t7 + t18 * t5)
      t19 = za(i2, i4)
      t20 = zb(i4, i2)
      t21 = t5 * t3
      t22 = t9 * (-2 * t13 * t21 * t8 - t2)
      t12 = 0.1D1 / t12
      t7 = 0.1D1 / t7
      t23 = 0.1D1 / t16
      t24 = 0.1D1 / t9
      t18 = 0.1D1 / t18
      t14 = 0.1D1 / t14
      t25 = 0.1D1 / t5
      t26 = 0.1D1 / t6
      t27 = t2 * t17
      t28 = t19 * t20
      t29 = t15 + t27 + t28
      t30 = t4 ** 2
      t31 = t14 ** 2
      t32 = t14 * t31
      t33 = t3 ** 2
      t34 = t3 * t33
      t35 = mt ** 2
      t36 = t1 ** 2
      t37 = t18 ** 2
      t38 = t5 * t7
      t39 = t22 * t23
      t40 = t3 * t26
      t41 = t5 * t33
      t42 = t33 * t26
      t43 = t12 * t8
      t44 = t7 * t8
      t45 = t24 * (t28 + t27) + t8
      t46 = t16 * t24
      t47 = t46 + t8
      t48 = t19 * t16
      t49 = t40 + t38
      t50 = t8 ** 2
      t51 = t8 * t50
      t52 = t2 ** 2
      t53 = t1 * t19 * t7
      t54 = t14 * t12
      t55 = t12 * t16
      t56 = t7 * t26
      t57 = t56 * t54
      t58 = t46 * t4
      t59 = t3 * t14
      t60 = t40 * t2
      t17 = t11 * (-t56 * t11 * t5 * t51 * t34 * t9 * t16 * t32 + t54 * 
     #t33 * (t56 * t36 * t2 * t30 * t9 * t14 + t21 * (t49 * t14 * t16 - 
     #t56) + t11 * t14 * (t56 * t21 * t16 + t49 * t9 * t2)) * t50 + t1 *
     # (t33 * (-t57 * t2 * t4 + t55 * t2 * t4 * t7 * (t11 * t26 + t5) * 
     #t31) + t34 * (-t58 * t56 * t5 * (t19 ** 2 * t20 ** 2 + t17 ** 2 * 
     #t52) * t32 + t55 * (t38 * t19 + (t53 + t2) * t26 * t4) * t31 - t57
     # * t19) + t48 * t33 ** 2 * t12 * t31 * t26 + t58 * t6 * t7 * t18 *
     # t37) * t8 + t60 * t7 * t24 * (-t18 * t25 + t59) * t35)
      t19 = t13 ** 2
      t20 = t5 ** 2
      t54 = t12 ** 2
      t57 = t26 * t34
      t58 = t15 * t6 * t16
      t61 = t18 * t7
      t6 = t50 * t30 * t36 * (t20 * (t15 * t54 * t34 * t7 * t31 + t55 * 
     #t7 * t34 * t29 * t32) + t5 * (-t56 * t35 * t34 * t19 * t23 * t14 +
     # t57 * (t11 * t35 * t7 * t19 + t15 * (t11 * t7 + t3) * t54) * t31 
     #+ t57 * t55 * (t3 * (t15 * (-t10 * t12 + 1) + t27 + t28) + t11 * t
     #29 * t7 - t15 * t36 * t30 * t12 * t7) * t32) + t61 * (t54 * (-t15 
     #* t3 * t18 + t58 * t37) + t55 * t6 * t37 + t42 * t35 * t19 * t23) 
     #- t58 * t54 * t5 * t20 * t34 * t7 * t32)
      ret = 32 * t44 * t30 * t36 * (t3 * (t43 * t37 + (t43 * t29 * t26 *
     # t31 + t16 * (t28 * t8 + t27 * (t24 * t28 + t8)) * t26 * t32) * t5
     # * t33) + t24 * (-t22 * t37 * t25 + t41 * (t16 * t40 + t22) * t31 
     #+ t40 * (t25 * t39 - t7) * t18 + t42 * (-t39 + t38) * t14) * t13 *
     # t35) + 8 * t42 * t11 * t14 * (t8 * t2 * t3 * t13 + (t1 * (t48 * t
     #45 * t3 + t2 * (t8 * (t15 + t16) + t28 * t47 + t27 * t47) * t4) + 
     #t21 * t8 * t16 * t45) * t14 * t7) - 16 * t17 - 64 * t6 + 128 * t54
     # * t32 * t16 * t9 * t51 * t5 * t30 * t34 * t36 * (t11 * t49 + t21)
     # - 4 * t11 * (t59 * (t52 * t7 + t60 * (t25 * (-t53 + t2) - t46 * t
     #7) + t41 * t44 * t24 * t26) - t61 * (t42 * t8 * t24 + t25 * t52))

      hjetmass_triangle_ppmm_s23_s123_0_dp = ret/32d0/(0,1d0)
      return

      end function

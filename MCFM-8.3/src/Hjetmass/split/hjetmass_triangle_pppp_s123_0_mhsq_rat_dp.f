
      double complex function hjetmass_triangle_pppp_s123_0_mhsq_rat_dp 
     &     (i1,i2,i3,i4,za,zb)
      implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision cg

      t1 = za(i1, i4)
      t2 = za(i2, i4)
      t3 = zb(i4, i1)
      t4 = za(i1, i2)
      t5 = zb(i2, i1)
      t6 = za(i1, i3)
      t7 = zb(i3, i1)
      t8 = za(i2, i3)
      t9 = zb(i3, i2)
      t10 = t4 * t5
      t11 = t6 * t7
      t12 = t8 * t9
      t13 = t10 + t11 + t12
      t14 = zb(i4, i2)
      t15 = za(i3, i4)
      t16 = zb(i4, i3)
      t17 = t1 * t3
      t18 = t2 * t14
      t19 = t15 * t16
      t20 = t17 + t18 + t19
      if ( dreal(t20) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t20 = cg * cdsqrt(t20 ** 2) + t17 + t18 + t19
      t20 = 0.1D1 / t20
      t10 = -2 * t17 * t13 * t20 + t10 + t11
      t11 = t14 * t8 + t3 * t6
      t21 = -2 * t1 * t13
      t22 = t21 * t14 * t20 + t6 * t9
      t21 = t21 * t16 * t20 - t4 * t9
      t23 = t2 * t22
      t24 = t3 * (t15 * t21 + t23)
      t25 = -t16 * t8 + t3 * t4
      t26 = t14 * t4 + t16 * t6
      t27 = -2 * t15 * t13 * t3 * t20 - t5 * t8
      t28 = t1 * (t14 * (-2 * t2 * t13 * t3 * t20 + t7 * t8) + t16 * t27
     #)
      t13 = 0.1D1 / t13
      t29 = 0.1D1 / t6
      t30 = 0.1D1 / t4
      t24 = 0.1D1 / t24
      t27 = 0.1D1 / t27
      t8 = 0.1D1 / t8
      t31 = 0.1D1 / t1
      t28 = 0.1D1 / t28
      t32 = 0.1D1 / t21
      t15 = 0.1D1 / t15
      t33 = t28 + t24
      t34 = t2 * t26
      t35 = t1 * t25
      t33 = -t34 * t33 + t35 * t33
      t36 = t10 * t14
      t37 = t36 * t2
      t27 = t10 * t27
      t38 = t16 * t30
      t39 = t38 * t8
      t40 = t1 * t10
      t41 = t14 * t29
      t42 = t30 * t29
      t43 = t14 * t25
      t44 = t16 * t11
      t45 = t44 + t43
      t46 = t5 * t29
      t47 = t7 * t30
      t48 = t47 + t46
      t49 = t3 ** 2
      t50 = t30 ** 2
      t51 = t30 * t1 * t15
      t21 = t21 * t30 * t8
      t52 = t13 * t26
      t2 = t20 * (t2 * (t49 * (t30 * t26 * ((-t51 - t8) * t29 * t22 - t2
     #1) * t24 + t40 * t30 * (-t21 * t45 + (-t45 * t8 + t51 * (-t44 - t4
     #3)) * t29 * t22) * t24 ** 2) - t14 * t26 * t13 * t8 * t31 * t48 - 
     #t17 * t16 * t10 * t22 * t29 * t50 * t32 * t15 * (t1 * t11 * t15 + 
     #t26) * t24) - t19 * t26 * t13 * t8 * t31 * t48 + t52 * (t8 * t48 *
     # t3 + t38 * (t30 * (-t12 * t29 - t7) - t46)))
      t12 = (t8 + t51) * t29
      t17 = t27 * t14
      t43 = t42 * (t15 * (t9 * t26 * t32 - t14) - t17 * t31)
      t44 = 16
      t45 = 32
      ret = t45 * ((t8 * (t5 * (t41 * t4 + t16) + t7 * (t38 * t6 + t14))
     # - t42 * (t18 + t19) * t20 * t26 * t9) * t31 * t13 + t20 * t30 * (
     #t3 * (t29 * (t30 * (t26 * ((-t37 * t24 + 1) * t15 * t1 + t27 * (-t
     #37 * t28 + 1)) + t35 * t14 * t10 * (t1 * t24 * t15 + t27 * t28)) +
     # t36 * t33 * t8) + t39 * t10 * t33) - t41 * t34 * t30 * t15 + t42 
     #* t15 * t16 * t1 * (t10 * t26 + (t40 - t23) * t15 * t11) * t32)) +
     # 64 * t2 + t44 * (t3 * (t28 * t10 * (t39 + t41 * (t27 * t30 + t8))
     # + (-t12 * t5 - t47 * t8) * t24 * t26) + t31 * (t41 + t38) * t13 *
     # t9 + t50 * t29 * ((-t34 + t35) * t15 * t32 * t22 * t16 + t17 * (-
     #t34 * t31 + t25)) * t20 + t24 * (t12 * t22 + t21) * t49) - 8 * t43
     # + 96 * t52 * t42 * t9 * t3 * t20

      hjetmass_triangle_pppp_s123_0_mhsq_rat_dp = ret/32d0/(0,1d0)
      return

      end function


      double complex function hjetmass_triangle_pppm_s134_0_mhsq_rat_dp 
     &     (i1,i2,i3,i4,za,zb)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision cg

      t1 = zb(i2, i1)
      t2 = zb(i3, i1)
      t3 = za(i1, i4)
      t4 = za(i3, i4)
      t5 = zb(i3, i2)
      t6 = t1 * t3 - t4 * t5
      t7 = za(i1, i2)
      t8 = za(i2, i3)
      t9 = za(i2, i4)
      t10 = zb(i4, i2)
      t11 = t7 * t1
      t12 = t8 * t5
      t13 = t9 * t10
      t14 = t12 + t13 + t11
      if ( dreal(t14) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t12 = cg * cdsqrt(t14 ** 2) + t11 + t12 + t13
      t14 = zb(i4, i1)
      t15 = zb(i4, i3)
      t16 = za(i1, i3)
      t17 = t16 * t2
      t18 = t3 * t14
      t19 = t4 * t15
      t20 = t19 + t17 + t18
      t12 = 0.1D1 / t12
      t21 = 2 * t9
      t22 = t21 * t1 * t20 * t12 - t2 * t4
      t11 = -2 * t11 * t20 * t12 + t17 + t18
      t18 = t10 * t22
      t23 = t1 * t11
      t24 = t7 * (t23 - t18)
      t25 = 2 * t8
      t26 = t25 * t1 * t20 * t12 + t14 * t4
      t18 = t7 * (t26 * t5 + t18)
      t27 = t10 * t3 + t16 * t5
      t28 = 2 * t7
      t29 = t28 * t5 * t20 * t12 + t15 * t3
      t28 = t9 * (t28 * t10 * t20 * t12 - t15 * t16)
      t30 = t1 * (t11 * t7 - t28)
      t21 = -t21 * t5 * t20 * t12 + t2 * t3
      t28 = t1 * (t29 * t8 + t28)
      t16 = t1 * t16 + t10 * t4
      t8 = 0.1D1 / t8
      t28 = 0.1D1 / t28
      t4 = 0.1D1 / t4
      t3 = 0.1D1 / t3
      t15 = 0.1D1 / t15
      t30 = 0.1D1 / t30
      t31 = 0.1D1 / t7
      t14 = 0.1D1 / t14
      t32 = t9 ** 2
      t33 = t1 * t4
      t24 = 0.1D1 / t24
      t34 = 0.1D1 / t29
      t26 = 0.1D1 / t26
      t18 = 0.1D1 / t18
      t35 = 0.1D1 / t11
      t36 = t28 + t30
      t37 = t1 ** 2
      t38 = t30 ** 2
      t39 = t28 ** 2
      t40 = t11 ** 2
      t41 = t9 * t29 * t31
      t42 = t2 * t5 * t6
      t43 = 0.1D1 / t20
      ret = 64 * t33 * t32 * t12 * t3 * (t5 * t27 * t21 * t31 * t15 * t3
     #0 - t23 * t6 * t8 * t14 * t28) + 32 * t3 * t4 * t12 * (t7 * (t42 *
     # t1 * t15 * t14 * (t24 + t18) * t22 + t6 * (-t5 ** 2 * t24 * t15 *
     # t35 + t37 * t26 * t18 * t14) * t22 ** 2) + t9 * t1 * (t1 * (-t32 
     #* t40 * t8 * t39 * t14 * (t10 * t6 + t16 * t5) + (t41 * t30 * (t35
     # * (-t41 - t21) + t13 * t21 * t30) + (t21 * t36 + ((t39 - t38) * t
     #21 * t10 + t5 * (-t25 * t5 * t20 * t12 + t17 + t19) * t39) * t11 *
     # t9) * t14 * t2) * t15 * t27) - t42 * t11 * t14 * t15 * t36 + t9 *
     # t27 * t29 * t15 * t38 * (t2 * t11 * t14 - t41) * t37 + t5 * t9 * 
     #t40 * t8 * t28 * t14 * (-t9 * t16 * t8 + t6) * t34)) - 128 * t2 **
     # 2 * t6 * t12 * t14 * t15 * t43 * (-t3 * t5 + t33)

      hjetmass_triangle_pppm_s134_0_mhsq_rat_dp = ret/32d0/(0,1d0)
      return

      end function

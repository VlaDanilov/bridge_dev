
      double complex function hjetmass_qqbgg_triangle_pmmp_s12_s124_0_dp
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

      t1 = za(i1, i2)
      t2 = zb(i4, i1)
      t3 = zb(i4, i2)
      t4 = zb(i2, i1)
      t5 = zb(i4, i3)
      t6 = zb(i3, i2)
      t7 = za(i1, i4)
      t8 = za(i2, i4)
      t9 = t7 * t2
      t10 = t8 * t3
      t11 = t10 + t9
      if ( dreal(t11) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t11 = cg * cdsqrt(t11 ** 2) + t10 + t9
      t12 = 0.1D1 / t11
      t13 = -2 * t7
      t14 = t13 * t4
      t15 = t1 * (t14 * t5 * t12 - t6)
      t16 = za(i2, i3)
      t17 = t1 * t4
      t13 = t17 * (t13 * t2 * t12 + 1)
      t18 = za(i3, i4)
      t14 = t10 * t14 * t1 * t2 * t12 + t2 * (t13 * t7 + t15 * t18)
      t19 = za(i1, i3)
      t20 = zb(i3, i1)
      t21 = t16 * t3
      t22 = t19 * t2
      t15 = 0.1D1 / t15
      t14 = 0.1D1 / t14
      t23 = 0.1D1 / t18
      t24 = 0.1D1 / t4
      t25 = 0.1D1 / t5
      t26 = t2 ** 2
      t27 = t2 * t26
      t28 = t16 ** 2
      t29 = t7 * t16
      t30 = t8 * t19
      t31 = t13 * t24
      t32 = t1 * t26
      t33 = t2 * t14
      t11 = 0.1D1 / t11
      t34 = t6 * t16
      t35 = t19 * t20
      t36 = t35 + t34
      t37 = t14 ** 2
      t38 = t14 * t37
      t39 = t1 ** 2
      t40 = t15 ** 2
      t41 = mt ** 2
      t42 = t18 * t13
      t43 = t23 * t15
      t44 = t26 * t14
      t45 = t1 * t39 * t4
      t46 = t9 * t25 + t18
      t47 = t10 * t25
      t48 = t47 + t18
      t49 = t3 ** 2
      t50 = t8 ** 2
      t51 = t11 * t13
      t52 = t21 * t11
      t4 = t3 * (t8 * (t26 * (t52 * t13 * t48 * t37 - t52 * t25 * t14) +
     # t27 * (-t19 * t11 * t25 * t14 + t51 * (t18 * t19 + (t30 + t29) * 
     #t25 * t3) * t37 - t42 * t3 * t24 * t25 * (t19 ** 2 * t20 ** 2 + t2
     #8 * t6 ** 2) * t38) + t31 * t3 * t5 * t15 * t40 + t51 * t26 ** 2 *
     # t7 * t19 * t37 * t25) * t1 + t44 * t11 * (t18 * t2 * (t46 * t14 *
     # t13 - t25) + t10 * t14 * (t46 * t16 * t4 + t42 * t2 * t25) + t49 
     #* t4 * t50 * t16 * t14 * t25) * t39 - t47 * t45 * t42 * t27 * t38 
     #+ t16 * t2 * t24 * t25 * (-t43 + t33) * t41)
      t6 = t35 + t17 + t34
      t7 = t18 ** 2
      t19 = t17 * t27 * t11 * t37
      t20 = t17 * t11
      t5 = t11 * t8 * t49 * t39 * (t18 * (t19 * t25 * (t10 + t9) + t25 *
     # t13 * t27 * (-t17 * t49 * t50 * t11 + t10 * t6 + t9 * (t17 * (-t9
     # * t11 + 1) + t34 + t35)) * t38) + t40 * (-t20 * t2 + (t20 + 1) * 
     #t15 * t13 * t5) + t7 * (t13 * t27 * t6 * t38 + t19) - t51 * t17 * 
     #t18 * t7 * t27 * t5 * t38)
      t6 = t24 * t36 + t1
      t7 = t31 + t1
      ret = -4 * t3 * (t15 * (t32 * t24 * t25 + t23 * t28) + t33 * (-t32
     # * t18 * t24 * t25 - t28 + t16 * (t23 * (t30 - t29) + t31) * t25 *
     # t2)) - 32 * t10 * (t45 * t11 * t18 * t27 * t3 * t25 * t37 + t24 *
     # t12 * (t21 + t22) * t41 * (t13 * t23 * t40 - t43 * t2 * t25 + t44
     # * (-t42 * t14 + t25)) + t2 * (t11 * t40 + (t11 * t36 * t25 * t37 
     #+ t13 * t36 * t25 * t38) * t18 * t26) * t3 * t39 + t34 * t35 * t42
     # * t1 * t27 * t3 * t24 * t25 * t38) + 16 * t4 + 64 * t5 - 8 * t44 
     #* t25 * t3 * (t9 * t1 * t16 * t12 + t14 * (t22 * t13 * t6 + t21 * 
     #(t1 * (t17 + t13) + t35 * t7 + t34 * t7)) * t8 + t42 * t33 * t1 * 
     #t6) - 128 * t45 * t42 * t38 * t11 ** 2 * t8 * t49 * t27 * (t10 * t
     #18 + t9 * t48)

      hjetmass_qqbgg_triangle_pmmp_s12_s124_0_dp = ret/32d0/(0,1d0)
      return

      end function

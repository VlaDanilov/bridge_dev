
      double complex function hjetmass_triangle_ppmm_s134_0_mhsq_rat_dp 
     &     (i1,i2,i3,i4,za,zb)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision cg

      t1 = za(i2, i4)
      t2 = zb(i2, i1)
      t3 = za(i1, i4)
      t4 = za(i3, i4)
      t5 = zb(i3, i2)
      t6 = t2 * t3 - t4 * t5
      t7 = za(i1, i2)
      t8 = za(i2, i3)
      t9 = zb(i4, i2)
      t10 = t7 * t2
      t11 = t8 * t5
      t12 = t1 * t9
      t13 = t11 + t12 + t10
      if ( dreal(t13) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t11 = cg * cdsqrt(t13 ** 2) + t10 + t11 + t12
      t12 = zb(i3, i1)
      t13 = zb(i4, i3)
      t14 = zb(i4, i1)
      t15 = za(i1, i3)
      t16 = t15 * t5 + t3 * t9
      t17 = t2 * t15 + t4 * t9
      t18 = t3 * t14
      t19 = t15 * t12
      t20 = t13 * t4 + t18 + t19
      t11 = 0.1D1 / t11
      t21 = 2 * t8 * t2 * t20 * t11 + t14 * t4
      t18 = -2 * t10 * t20 * t11 + t18 + t19
      t19 = 2 * t1 * t2 * t20 * t11 - t12 * t4
      t22 = t9 * t19
      t23 = t7 * (t21 * t5 + t22)
      t24 = t2 * t18
      t22 = t7 * (t24 - t22)
      t25 = 2 * t7
      t15 = t1 * (t25 * t9 * t20 * t11 - t13 * t15)
      t26 = t2 * (t18 * t7 - t15)
      t25 = t25 * t5 * t20 * t11 + t13 * t3
      t15 = t2 * (t25 * t8 + t15)
      t20 = 0.1D1 / t20
      t15 = 0.1D1 / t15
      t14 = 0.1D1 / t14
      t3 = 0.1D1 / t3
      t13 = 0.1D1 / t13
      t27 = 0.1D1 / t12
      t28 = 0.1D1 / t7
      t29 = 0.1D1 / t1
      t30 = t24 * t15
      t31 = t30 + t28
      t32 = t1 * t6 * t11 * t27
      t33 = t13 * t3 * t2
      t25 = 0.1D1 / t25
      t34 = 0.1D1 / t8
      t26 = 0.1D1 / t26
      t23 = 0.1D1 / t23
      t35 = 0.1D1 / t18
      t36 = t10 * t23 + t35
      t37 = t18 ** 2
      t38 = t15 ** 2
      t39 = t2 ** 2
      t40 = t26 ** 2
      t41 = t1 ** 2
      t42 = t16 * t28 ** 2
      t9 = t6 * t9
      t43 = t8 * t14
      t44 = t13 * t15
      t45 = t44 * t2
      t46 = t15 * t37 * t1
      t47 = t3 * t11 * t2
      t48 = 0.1D1 / t5
      t49 = 0.1D1 / t21
      t22 = 0.1D1 / t22
      t50 = t8 * t19 * t29
      t21 = t50 + t21
      t22 = t22 * t35
      t51 = t18 * t25
      t52 = t26 * t28
      t7 = t47 * t14 * (t6 * (t10 * (-t22 * t19 * t21 * t27 - t23 * t13 
     #* t21) - t19 * t48 * t23 * t27 * (t50 * t49 + 1) * t7 * t39 - t35 
     #* t13 * t21) + t13 * t31 * t17 * t1 - t24 * t27 * (t51 * t34 * t15
     # + t52) * t17 * t41)
      t3 = t3 * t14
      t21 = t33 * t4 * t6 * t11 * t20 * (t43 * t12 * t29 + 1)
      ret = -64 * t33 * (-t4 ** 2 * t12 * t29 * t14 * t20 + t32 * t31) -
     # 48 * t33 * t8 * t6 * t11 * t14 * (t30 + t28) - 32 * t47 * (t46 * 
     #(t45 * t8 * t14 + (t18 * (-t2 * t15 * t14 * t25 - t34 * t14 * t25 
     #** 2) + t45) * t27 * t1) * t17 * t5 + t27 * (t2 * (-t25 * t9 * t38
     # * t14 * t18 * t37 + t9 * (t8 * t40 * t14 * t28 + t13 * t38) * t37
     # - t43 * t42 * t26 * t18) - t39 * t8 * t16 * t37 * t40 * t14 * t28
     # + t42 * t13) * t41 - t6 * t13 * t27 * t36 * t19 + t43 * t13 * (t9
     # * t2 * t37 * t38 + t42) * t1) - 16 * t7 + 8 * t3 * t4 * t2 * (-t1
     #3 * t28 + t2 * (t18 * (t52 * t1 * t27 - t44) + t46 * t25 * t34 * t
     #27) + t13 * t36 * t29 * t19 + t10 * t27 * (t2 * t48 * t49 * t23 + 
     #t22) * t29 * t19 ** 2) + 80 * t3 * t32 * t39 * t18 * (t51 * t15 + 
     #t52 * t8) + 128 * t21

      hjetmass_triangle_ppmm_s134_0_mhsq_rat_dp = ret/32d0/(0,1d0)
      return

      end function


      double complex function hjetmass_box_ppmm_0_0_s12_mhsq_s34_s124_dp 
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision mt

      t1 = za(i1, i2)
      t2 = za(i3, i4)
      t3 = zb(i2, i1)
      t4 = zb(i4, i3)
      t5 = za(i1, i4)
      t6 = zb(i4, i1)
      t7 = za(i2, i4)
      t8 = zb(i4, i2)
      t9 = zb(i3, i1)
      t10 = zb(i3, i2)
      t11 = t5 * t9
      t12 = t10 * t7 + t11
      t13 = mt ** 2
      t14 = za(i1, i3)
      t15 = za(i2, i3)
      t16 = t14 * t6
      t17 = t15 * t8
      t18 = t17 + t16
      t19 = t2 * t4
      t20 = t19 * t18
      t21 = t1 * t3
      t22 = t5 * t6
      t23 = t7 * t8
      t24 = t19 * (t23 + t21 + t22)
      t25 = 16 * t13 * t12 * t20 + 4 * t24 ** 2
      t25 = cdsqrt(t25)
      t24 = 2 * t24
      t26 = -t24 - t25
      t20 = 0.1D1 / t20
      t16 = t17 + t16
      t16 = t19 * t20 * (t21 * t16 + t22 * t16 + t23 * t16)
      t17 = (0.1D1 / 0.2D1)
      t18 = t17 * t25 * t20 * t18
      t19 = 2 * t23 + 2 * t21 + 2 * t22
      t21 = t19 - t18 - t16
      t12 = 0.1D1 / t12
      t11 = t11 * t17
      t22 = t6 * (t14 * t26 * t20 / 4 + t5) - t11 * t21 * t12
      t24 = -t24 + t25
      t16 = t19 + t18 - t16
      t11 = t6 * (t14 * t24 * t20 / 4 + t5) - t11 * t16 * t12
      t18 = 0.1D1 / t22
      t4 = 0.1D1 / t4
      t5 = 0.1D1 / t5
      t19 = 0.1D1 / t6
      t1 = 0.1D1 / t1
      t7 = 0.1D1 / t7
      t11 = 0.1D1 / t11
      t22 = t8 * t9
      t25 = t22 * t19
      t27 = t6 * t10
      t28 = -t22 + t27
      t9 = t2 * t9 * t7 + t28 * t4
      t29 = t24 ** 2
      t30 = t26 ** 2
      t31 = t30 * t18
      t32 = t29 * t11
      t33 = t32 + t31
      t34 = t2 ** 2
      t35 = t2 * t34
      t10 = t1 * (t4 * (t23 * (-t25 + t10) * t5 + t27 - t22) + t5 * (t25
     # - t10) * t2)
      t25 = t3 * t11
      t36 = t3 * t18
      t37 = t2 * t6
      t38 = t20 * t34
      t39 = t8 * t4
      t22 = t28 * t7 * t2 + t39 * (t22 - t27)
      t27 = t20 ** 2
      t28 = t12 * t5
      t40 = t11 * t24 + t18 * t26
      t41 = t26 + t24
      t42 = t41 * t5
      t43 = t41 * t1
      t44 = t30 * t21
      t45 = t29 * t16
      ret = t17 * t38 * (t12 * (t16 * t24 * (t25 * t9 + t10) + (t36 * t9
     # + t10) * t21 * t26) + t37 * t8 * t4 * t5 * t7 * t33 * t20 * t13) 
     #+ t28 * t27 * t35 * (t32 * t16 * t22 + t31 * t22 * t21) / 8 + t28 
     #* t20 * t27 * t34 ** 2 * t6 * t8 * t7 * (t24 * t29 * t16 * t11 + t
     #26 * t30 * t21 * t18) / 16 - 8 * t1 * t3 * t2 * t13 * ((-t23 * t4 
     #+ t2) * t19 * t5 - t4) + t38 * (t3 * (t1 * (-t15 * t7 * t41 + t42 
     #* t14) + t5 * t7 * t40 * t2 * t13 - t39 * t13 * t5 * t40) + t42 * 
     #t2 * t8 * t1 + t43 * t37 * t7) - t27 * t35 * (t28 * t8 * t1 * (t45
     # + t44) + (t12 * (t45 * (t25 + t1) + t44 * (t36 + t1)) + t3 * t33 
     #* t5 * t14 + t5 * t33 * t8 * t2) * t7 * t6) / 4 - 2 * t38 * t4 * t
     #13 * (t42 * t8 * t1 + (t3 * t40 + t43) * t7 * t6)

      hjetmass_box_ppmm_0_0_s12_mhsq_s34_s124_dp = ret/32d0/(0,1d0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_box_pmmp_0_0_s12_mhsq_s34_s123
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt

      t1 = za(i1, i2)
      t2 = zb(i2, i1)
      t3 = za(i1, i3)
      t4 = zb(i3, i1)
      t5 = za(i2, i3)
      t6 = zb(i3, i2)
      t7 = za(i3, i4)
      t8 = zb(i4, i1)
      t9 = zb(i4, i3)
      t10 = t3 * t8
      t11 = t5 * zb(i4, i2)
      t12 = t11 + t10
      t13 = t7 * t9
      t14 = t13 * t12
      t15 = t1 * t2
      t3 = t3 * t4
      t16 = t5 * t6
      t17 = t13 * (t16 + t15 + t3)
      t18 = za(i2, i4)
      t6 = t4 * za(i1, i4) + t6 * t18
      t19 = mt ** 2
      t20 = 16 * t19 * t6 * t14 + 4 * t17 ** 2
      t20 = sqrt(t20)
      t14 = 0.1q1 / t14
      t21 = t16 + t3
      t22 = t13 * t14
      t10 = t22 * (t10 * t21 + t11 * t21 + t15 * (t11 + t10))
      t11 = 2
      t3 = t11 * (t16 + t15 + t3)
      t12 = t20 * t14 * t12 / 2
      t14 = t3 - t10 - t12
      t3 = t3 - t10 + t12
      t10 = t11 * t17
      t2 = 0.1q1 / t2
      t1 = 0.1q1 / t1
      t6 = 0.1q1 / t6
      t12 = t14 + t3
      t15 = t14 ** 2
      t16 = t3 ** 2
      t17 = t6 ** 2
      ret = t11 * t6 * (t5 ** 2 * t9 * t1 * t12 + (t12 * t8 ** 2 + (-t14
     # * t15 - t3 * t16) * t17 * t1 * t18 * t9 * t4) * t2 * t7) + 4 * t6
     # * t1 * t2 * t5 * (t19 * t8 * t12 + t13 * (t16 + t15) * t6 * t4) +
     # t22 * t17 * t8 * t5 * t1 * t2 * (t15 * (-t10 - t20) + t16 * (-t10
     # + t20)) / 2

      hjetmass_qqbgg_box_pmmp_0_0_s12_mhsq_s34_s123 = ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_box_pmmp_0_0_s12_mhsq_s34_s124
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt

      t1 = za(i1, i2)
      t2 = za(i3, i4)
      t3 = zb(i2, i1)
      t4 = zb(i4, i3)
      t5 = za(i1, i4)
      t6 = zb(i4, i1)
      t7 = za(i2, i4)
      t8 = zb(i4, i2)
      t9 = za(i2, i3)
      t10 = t6 * za(i1, i3) + t8 * t9
      t11 = mt ** 2
      t12 = zb(i3, i1)
      t13 = t5 * t12
      t14 = t7 * zb(i3, i2)
      t15 = t14 + t13
      t16 = t2 * t4
      t17 = t16 * t15
      t18 = t1 * t3
      t5 = t5 * t6
      t8 = t7 * t8
      t19 = t16 * (t8 + t18 + t5)
      t20 = 16 * t11 * t10 * t17 + 4 * t19 ** 2
      t20 = sqrt(t20)
      t19 = 2 * t19
      t21 = -t20 - t19
      t19 = t20 - t19
      t17 = 0.1q1 / t17
      t22 = t8 + t5
      t13 = t16 * t17 * (t13 * t22 + t18 * (t14 + t13) + t14 * t22)
      t14 = t20 * t17 * t15 / 2
      t5 = 2 * t8 + 2 * t18 + 2 * t5
      t1 = 0.1q1 / t1
      t8 = 0.1q1 / t10
      t3 = 0.1q1 / t3
      t10 = t21 ** 2
      t15 = t19 ** 2
      t18 = t21 + t19
      ret = t16 * t3 * t17 ** 2 * t1 * ((t10 * (-t13 + t5 - t14) + t15 *
     # (-t13 + t5 + t14)) * t8 * t9 * t6 - t17 * t12 * t7 * (t21 * t10 +
     # t19 * t15)) / 4 - t17 * (-t2 * t6 ** 2 * t3 * t18 + (-t18 * t9 **
     # 2 + (t15 + t10) * t3 * t17 * t7 * t6 * t2) * t1 * t4) + 2 * t11 *
     # t17 * t9 * t6 * t1 * t3 * (t21 + t19)

      hjetmass_qqbgg_box_pmmp_0_0_s12_mhsq_s34_s124 = ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_box_pmmp_0_s12_0_mhsq_s123_s124
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt

      t1 = za(i1, i3)
      t2 = za(i1, i4)
      t3 = zb(i3, i1)
      t4 = zb(i4, i1)
      t5 = za(i2, i3)
      t6 = zb(i3, i2)
      t7 = za(i2, i4)
      t8 = zb(i4, i2)
      t9 = za(i1, i2)
      t10 = za(i3, i4)
      t11 = zb(i2, i1)
      t12 = zb(i4, i3)
      t13 = t1 * t4
      t14 = t5 * t8
      t15 = t14 + t13
      t16 = t10 ** 2 * t12 ** 2
      t17 = t16 * t15
      t18 = t2 * t4
      t19 = t7 * t8
      t20 = t10 * t12
      t21 = t20 * t9 * t11
      t22 = t5 * t6
      t23 = t20 * ((-t19 - t18) * t3 * t1 + t21 - t22 * (t19 + t18))
      t24 = t2 * t3 + t6 * t7
      t25 = mt ** 2
      t23 = 16 * t20 * t25 * t24 * t17 + 4 * t23 ** 2
      t17 = 0.1q1 / t17
      t26 = t14 + t13
      t27 = t1 * t26 * t3 + t22 * t26
      t28 = t19 + t18
      t22 = t28 * t3 * t1 + t22 * t28 - t21
      t23 = t20 * sqrt(t23) * t17 * t15 / 4
      t16 = t16 * t17 * (t18 * t27 + t19 * t27 + t21 * (-t14 - t13)) / 2
      t17 = -t22 - t23 + t16
      t21 = t20 * t24
      t16 = -t22 + t23 + t16
      t21 = 0.1q1 / t21
      t22 = 0.1q1 / t10
      t18 = (t20 + t18 + t19) * t22
      t19 = t18 * t5
      t20 = t17 * t21
      t23 = t20 * t7 * t12 + t19
      t24 = t16 * t21
      t7 = t24 * t7 * t12 + t19
      t18 = t18 * t1
      t19 = 0.1q1 / t12
      t11 = 0.1q1 / t11
      t9 = 0.1q1 / t9
      t27 = t17 + t16
      t28 = t16 ** 2
      t29 = t17 ** 2
      t30 = t22 * t19
      t31 = t27 * t21
      t32 = t5 ** 2
      t33 = -t16 * t7 - t17 * t23
      t16 = t17 + t16
      t17 = t16 * t5 + t19 * t33
      t34 = 8
      ret = t34 * (t9 * t5 * (t5 * t22 * t15 + t21 * t33) + t21 * (t9 * 
     #(t13 * t17 + t14 * t17) - t4 * (t29 + t28) * t21 * t10) * t11 * t3
     #) + 12 * t21 * t32 * t12 * t9 * t16 + 16 * t9 * t11 * (t21 ** 2 * 
     #(-t23 * t29 - t28 * t7 + (t29 + t28) * t12 * t5) * t10 * t3 + t5 *
     # t4 * t25 * (t30 * t26 + t31)) - 4 * t4 * (-t31 * t10 * t4 * t11 -
     # t30 * t9 * (t20 * t2 * t12 + t24 * t2 * t12 + 2 * t18) * t32 + ((
     #t27 * t8 * t3 - t4 * t6 * t27) * t11 * t21 + t1 * t22 * t9 * (t23 
     #+ t7)) * t19 * t5)


      hjetmass_qqbgg_box_pmmp_0_s12_0_mhsq_s123_s124 = ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_box_pmpm_0_0_s12_mhsq_s34_s123
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt

      t1 = za(i1, i2)
      t2 = za(i3, i4)
      t3 = zb(i2, i1)
      t4 = zb(i4, i3)
      t5 = za(i1, i3)
      t6 = zb(i3, i1)
      t7 = za(i2, i3)
      t8 = zb(i3, i2)
      t9 = t1 * t3
      t10 = t5 * t6
      t11 = t7 * t8
      t12 = t2 * t4
      t13 = t12 * (t11 + t9 + t10)
      t14 = za(i2, i4)
      t8 = t6 * za(i1, i4) + t8 * t14
      t15 = mt ** 2
      t16 = zb(i4, i1)
      t5 = t5 * t16
      t17 = t7 * zb(i4, i2)
      t18 = t17 + t5
      t19 = t12 * t18
      t20 = 16 * t15 * t8 * t19 + 4 * t13 ** 2
      t20 = sqrt(t20)
      t13 = 2 * t13
      t21 = -t13 - t20
      t13 = -t13 + t20
      t19 = 0.1q1 / t19
      t5 = t17 + t5
      t5 = t12 * t19 * (t10 * t5 + t11 * t5 + t9 * t5)
      t9 = 2 * t11 + 2 * t9 + 2 * t10
      t10 = t20 * t19 * t18 / 2
      t3 = 0.1q1 / t3
      t1 = 0.1q1 / t1
      t11 = t21 + t13
      t17 = t21 ** 2
      t18 = t13 ** 2
      t8 = 0.1q1 / t8
      ret = t12 * t1 * t3 * t19 ** 2 * (-t6 * t8 * t14 * (t17 * (-t10 - 
     #t5 + t9) + t18 * (t10 - t5 + t9)) + (t13 * t18 + t21 * t17) * t19 
     #* t16 * t7) / 4 + t19 * (t2 * t6 ** 2 * t3 * t11 + (t11 * t14 ** 2
     # + (t18 + t17) * t3 * t19 * t7 * t6 * t2) * t1 * t4) - 2 * t15 * t
     #19 * t14 * t6 * t1 * t3 * (t21 + t13)

      hjetmass_qqbgg_box_pmpm_0_0_s12_mhsq_s34_s123 = ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_box_pmpm_0_0_s12_mhsq_s34_s124
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt

      t1 = za(i1, i2)
      t2 = zb(i2, i1)
      t3 = za(i1, i4)
      t4 = zb(i4, i1)
      t5 = za(i2, i4)
      t6 = zb(i4, i2)
      t7 = za(i3, i4)
      t8 = zb(i3, i1)
      t9 = zb(i4, i3)
      t10 = t3 * t8
      t11 = t5 * zb(i3, i2)
      t12 = t11 + t10
      t13 = t7 * t9
      t14 = t13 * t12
      t15 = za(i2, i3)
      t16 = t4 * za(i1, i3) + t6 * t15
      t17 = mt ** 2
      t18 = t1 * t2
      t3 = t3 * t4
      t6 = t5 * t6
      t19 = t13 * (t6 + t18 + t3)
      t20 = 16 * t17 * t16 * t14 + 4 * t19 ** 2
      t20 = sqrt(t20)
      t14 = 0.1q1 / t14
      t21 = t6 + t3
      t22 = t13 * t14
      t10 = t22 * (t10 * t21 + t11 * t21 + t18 * (t11 + t10))
      t11 = 2
      t12 = t20 * t14 * t12 / 2
      t3 = t11 * (t6 + t18 + t3)
      t6 = t3 - t12 - t10
      t3 = t3 + t12 - t10
      t10 = t11 * t19
      t1 = 0.1q1 / t1
      t2 = 0.1q1 / t2
      t12 = 0.1q1 / t16
      t14 = t3 ** 2
      t16 = t12 ** 2
      t18 = t6 ** 2
      t19 = t6 + t3
      ret = t11 * t12 * (t5 ** 2 * t9 * t1 * t19 + (t19 * t8 ** 2 + (t3 
     #* t14 + t6 * t18) * t16 * t1 * t15 * t9 * t4) * t2 * t7) - t22 * t
     #16 * t8 * t5 * t1 * t2 * (t18 * (-t10 - t20) + t14 * (-t10 + t20))
     # / 2 - 4 * t12 * t2 * t1 * t5 * (t17 * t8 * t19 + t13 * (t18 + t14
     #) * t12 * t4)

      hjetmass_qqbgg_box_pmpm_0_0_s12_mhsq_s34_s124 = ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_box_pmpm_0_s12_0_mhsq_s123_s124
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt

      t1 = za(i2, i4)
      t2 = zb(i3, i1)
      t3 = za(i1, i4)
      t4 = zb(i4, i1)
      t5 = zb(i4, i2)
      t6 = za(i3, i4)
      t7 = zb(i4, i3)
      t8 = t1 * t5
      t9 = t3 * t4
      t10 = t6 * t7
      t11 = zb(i3, i2)
      t12 = za(i1, i3)
      t13 = za(i2, i3)
      t14 = za(i1, i2)
      t15 = zb(i2, i1)
      t16 = t2 * t12
      t17 = t11 * t13
      t16 = t10 * (t10 * t14 * t15 - t8 * (t17 + t16) + t9 * (-t17 - t16
     #))
      t17 = t2 * t3
      t18 = t1 * t11
      t19 = t17 + t18
      t20 = mt ** 2
      t21 = t6 ** 2 * t7 ** 2 * (t12 * t4 + t13 * t5)
      t22 = 16 * t10 * t20 * t19 * t21 + 4 * t16 ** 2
      t22 = sqrt(t22)
      t16 = -2 * t16
      t23 = -t22 + t16
      t21 = 0.1q1 / t21
      t24 = 0.1q1 / t7
      t8 = (t9 + t10 + t8) * t24
      t9 = t8 * t11
      t10 = t23 * t21 * t6 / 4
      t16 = t22 + t16
      t22 = t16 * t21 * t6 / 4
      t8 = t8 * t2
      t25 = t10 * t4 - t8
      t4 = t22 * t4 - t8
      t8 = 0.1q1 / t14
      t14 = 0.1q1 / t15
      t15 = 0.1q1 / t6
      t26 = t16 * t4 + t23 * t25
      t27 = t23 + t16
      t28 = t15 * t26 + t2 * t27
      t29 = t14 * t21
      t30 = t23 ** 2
      t31 = t16 ** 2
      t32 = t1 * t7
      t33 = t8 * t21
      t34 = t2 ** 2
      ret = 16 * t20 * t1 * t2 * t8 * t15 * t14 * t24 * t19 + 3 * t29 * 
     #t6 * t34 * (t23 + t16) - 8 * t34 * t14 * t24 * t19 + t32 * t21 ** 
     #2 * t13 * t8 * (t31 + t30) / 2 - 2 * t29 * (-t2 * t26 + (t17 * t28
     # + t18 * t28) * t8 * t13) + t33 * (t1 * (t32 * t27 + (t27 * t12 * 
     #t1 - t3 * t13 * t27) * t15 * t2) + t29 * t13 * t7 * (t25 * t30 + t
     #31 * t4 + (t31 + t30) * t6 * t2)) - 4 * t14 * t2 * t1 * (t33 * t20
     # * t27 + (t11 * (t25 + t4) + t2 * (-t10 * t5 - t22 * t5 + 2 * t9))
     # * t15 * t24)

      hjetmass_qqbgg_box_pmpm_0_s12_0_mhsq_s123_s124 = ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_box_pmpp_0_0_s12_mhsq_s34_s123
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt

      t1 = zb(i3, i1)
      t2 = za(i1, i2)
      t3 = zb(i2, i1)
      t4 = za(i1, i3)
      t5 = za(i2, i3)
      t6 = zb(i3, i2)
      t7 = za(i3, i4)
      t8 = zb(i4, i1)
      t9 = zb(i4, i3)
      t10 = t4 * t8
      t11 = t5 * zb(i4, i2)
      t12 = t11 + t10
      t13 = t7 * t9
      t14 = t13 * t12
      t15 = t2 * t3
      t16 = t4 * t1
      t17 = t5 * t6
      t18 = t13 * (t17 + t15 + t16)
      t19 = za(i1, i4)
      t20 = za(i2, i4)
      t21 = t19 * t1
      t6 = t20 * t6 + t21
      t22 = mt ** 2
      t23 = 16 * t22 * t6 * t14 + 4 * t18 ** 2
      t23 = sqrt(t23)
      t14 = 0.1q1 / t14
      t24 = t17 + t16
      t11 = t13 * t14 * (t10 * t24 + t11 * t24 + t15 * (t11 + t10))
      t13 = 2
      t24 = (0.1q1 / 0.2q1)
      t12 = t24 * t23 * t14 * t12
      t15 = t13 * (t17 + t15 + t16)
      t17 = -t12 - t11 + t15
      t18 = t13 * t18
      t25 = -t23 - t18
      t6 = 0.1q1 / t6
      t26 = -(0.1q1 / 0.4q1)
      t10 = t10 * t26
      t27 = t21 * t24
      t28 = t10 * t25 * t14 + t27 * t17 * t6
      t11 = t12 - t11 + t15
      t12 = t23 - t18
      t10 = t10 * t12 * t14 + t27 * t11 * t6
      t15 = t19 * t5 - t20 * t4
      t18 = t12 * t14 * t11 * t6 * t9 * t15
      t15 = t25 * t14 * t17 * t6 * t9 * t15
      t4 = 0.1q1 / t4
      t23 = 0.1q1 / t2
      t27 = 0.1q1 / t19
      t3 = 0.1q1 / t3
      t29 = t4 * (t28 + t10)
      t30 = t3 * t9
      t31 = t8 * t3
      t32 = t1 * t3
      t33 = t12 * t11
      t34 = t17 + t11
      t35 = t4 * t1
      t36 = t30 * t6
      t37 = t35 * t36 * t23 * (t11 * t18 + t15 * t17)
      t38 = -t27 * (t32 * t7 + t20) + (-t31 * t19 * t23 + 1) * t4 * t5
      t3 = t9 * (t3 * t27 * t23 * (t1 * (t18 + t15) + t4 * (-t10 * t18 -
     # t15 * t28)) + (t25 * t17 * t38 + t33 * t38) * t14 * t6 * t9)
      t7 = t31 * t14 ** 2 * t9 ** 2 * t6 * t7 * t27 * (t11 * t12 ** 2 + 
     #t17 * t25 ** 2)
      ret = -t13 * t36 * ((t20 * (t1 * t34 + t4 * (-t10 * t11 - t17 * t2
     #8)) - t21 * t4 * t34 * t5) * t23 * t9 + t35 * (t11 * (-t22 + t10) 
     #+ t17 * (-t22 + t28) + t34 * t19 * t8)) - t24 * t3 + t26 * t37 - 8
     # * t9 * ((t2 * t27 - t31) * t4 * t22 + t32 * (-t9 * t5 * t23 + t8)
     #) + t7 / 8 + 4 * t30 * (t27 * (t1 * (-t28 - t10) + t4 * (t10 ** 2 
     #+ t10 * t22 + t22 * t28 + t28 ** 2)) + t29 * t8 - t29 * t23 * t9 *
     # t5) + t30 * (t9 * (t14 * ((t5 * (t28 * t4 - t1) + (-t16 + t28) * 
     #t27 * t20) * t25 * t17 + t33 * (t5 * (t10 * t4 - t1) + (-t16 + t10
     #) * t27 * t20)) * t23 * t6 - t21 * t20 * t4 * (t11 ** 2 + t17 ** 2
     #) * t23 * t6 ** 2) + t27 * t14 * t8 * (t12 * (t22 + t10) + t25 * (
     #t22 + t28) + t16 * (-t25 - t12)))

      hjetmass_qqbgg_box_pmpp_0_0_s12_mhsq_s34_s123 = ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_box_pmpp_0_0_s12_mhsq_s34_s124
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt

      t1 = za(i1, i2)
      t2 = zb(i2, i1)
      t3 = za(i1, i4)
      t4 = zb(i4, i1)
      t5 = za(i2, i4)
      t6 = zb(i4, i2)
      t7 = za(i3, i4)
      t8 = zb(i3, i1)
      t9 = zb(i4, i3)
      t10 = t3 * t8
      t11 = t5 * zb(i3, i2)
      t12 = t11 + t10
      t13 = t7 * t9
      t14 = t13 * t12
      t15 = za(i1, i3)
      t16 = za(i2, i3)
      t17 = t15 * t4
      t18 = t16 * t6 + t17
      t19 = mt ** 2
      t20 = t1 * t2
      t21 = t3 * t4
      t6 = t5 * t6
      t22 = t13 * (t6 + t20 + t21)
      t23 = 16 * t19 * t18 * t14 + 4 * t22 ** 2
      t23 = sqrt(t23)
      t14 = 0.1q1 / t14
      t24 = t6 + t21
      t11 = t13 * t14 * (t10 * t24 + t11 * t24 + t20 * (t11 + t10))
      t13 = 2
      t12 = t23 * t14 * t12 / 2
      t6 = t13 * (t6 + t20 + t21)
      t20 = -t12 + t6 - t11
      t22 = t13 * t22
      t24 = -t22 - t23
      t22 = -t22 + t23
      t6 = t12 + t6 - t11
      t11 = 0.1q1 / t18
      t10 = -t10 / 4
      t12 = t17 / 2
      t18 = t10 * t24 * t14 + t12 * t20 * t11
      t23 = -t15 * t5 + t16 * t3
      t25 = t20 * t24
      t26 = t25 * t14 * t11 * t9 * t23
      t27 = t22 * t6
      t23 = t27 * t14 * t11 * t9 * t23
      t10 = t10 * t22 * t14 + t12 * t6 * t11
      t12 = 0.1q1 / t1
      t28 = 0.1q1 / t15
      t2 = 0.1q1 / t2
      t3 = 0.1q1 / t3
      t29 = t27 + t25
      t30 = t8 * t2
      t31 = t2 * t28
      t32 = t20 + t6
      t17 = t17 * t3
      t33 = t3 * t4
      t34 = t2 * t11 * t9
      t35 = t2 * t9
      t36 = t31 * t14 ** 2 * t9 ** 2 * t11 * t7 * t8 * (t20 * t24 ** 2 +
     # t6 * t22 ** 2)
      t37 = t5 * t9 * t12
      t1 = t9 * ((-t1 * t28 + t30) * t3 * t19 - t4 * t2 * (t37 + t8))
      t38 = t3 * (t18 + t10)
      t37 = t35 * (t28 * (t3 * (t10 ** 2 + t10 * t19 + t18 ** 2 + t19 * 
     #t18) + t4 * (-t18 - t10)) + t38 * t8 + t38 * t37)
      ret = t13 * t34 * ((t16 * (t3 * (-t10 * t6 - t18 * t20) + t32 * t4
     #) - t17 * t32 * t5) * t12 * t9 + t33 * (t20 * (t19 - t18) + t6 * (
     #t19 - t10) + (-t20 - t6) * t15 * t8)) + t36 / 8 + 8 * t1 - t33 * t
     #34 * t12 * (t20 * t26 + t23 * t6) / 4 - t9 * ((t28 * (t16 * t29 + 
     #(-t27 - t25) * t2 * t7 * t4) + (t30 * t29 * t12 * t15 - t25 - t27)
     # * t3 * t5) * t14 * t11 * t9 + t31 * t12 * (t3 * (-t10 * t23 - t18
     # * t26) + t4 * (t26 + t23))) / 2 + t35 * (t9 * (t14 * (t25 * (t5 *
     # (-t18 * t3 + t4) + (t21 - t18) * t28 * t16) + t27 * (t5 * (-t10 *
     # t3 + t4) + (t21 - t10) * t28 * t16)) * t12 * t11 + t17 * t16 * (t
     #20 ** 2 + t6 ** 2) * t12 * t11 ** 2) + t28 * t14 * t8 * (t22 * (t1
     #9 + t10) + t24 * (t19 + t18) + t21 * (-t24 - t22))) + 4 * t37

      hjetmass_qqbgg_box_pmpp_0_0_s12_mhsq_s34_s124 = ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_box_pmpp_0_s12_0_mhsq_s123_s124
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt

      t1 = zb(i3, i1)
      t2 = zb(i4, i2)
      t3 = za(i1, i3)
      t4 = za(i1, i4)
      t5 = za(i3, i4)
      t6 = zb(i4, i1)
      t7 = zb(i4, i3)
      t8 = za(i2, i3)
      t9 = zb(i3, i2)
      t10 = za(i2, i4)
      t11 = za(i1, i2)
      t12 = zb(i2, i1)
      t13 = t4 * t6
      t14 = t10 * t2
      t15 = t5 * t7
      t16 = t15 * t11
      t17 = t16 * t12
      t18 = t8 * t9
      t19 = t15 * ((-t14 - t13) * t3 * t1 + t17 - t18 * (t14 + t13))
      t20 = t1 * t4 + t10 * t9
      t21 = mt ** 2
      t22 = t8 * t2
      t23 = t3 * t6
      t24 = t23 + t22
      t25 = t7 ** 2
      t26 = t5 ** 2 * t25
      t27 = t26 * t24
      t28 = 16 * t15 * t21 * t20 * t27 + 4 * t19 ** 2
      t28 = sqrt(t28)
      t19 = 2 * t19
      t29 = t28 - t19
      t27 = 0.1q1 / t27
      t30 = t23 + t22
      t30 = t3 * t30 * t1 + t18 * t30
      t31 = t14 + t13
      t18 = t31 * t3 * t1 + t18 * t31 - t17
      t24 = t15 * t28 * t27 * t24 / 4
      t17 = t26 * t27 * (t13 * t30 + t14 * t30 + t17 * (-t23 - t22)) / 2
      t22 = -t18 + t24 + t17
      t20 = t15 * t20
      t26 = 0.1q1 / t7
      t20 = 0.1q1 / t20
      t5 = 0.1q1 / t5
      t5 = (t15 + t13 + t14) * t5
      t13 = t5 * t26 * t3
      t14 = t22 * t20
      t15 = t14 * t4
      t26 = t15 + t13
      t30 = t23 / 4
      t31 = -t30 * t29 * t27 + t1 * t26
      t32 = t5 * t8
      t14 = t14 * t10 * t7 + t32
      t17 = -t18 - t24 + t17
      t18 = -t28 - t19
      t19 = t17 * t20
      t24 = t19 * t4
      t13 = t24 + t13
      t28 = -t18 * t30 * t27 + t1 * t13
      t10 = t19 * t10 * t7 + t32
      t5 = t5 * t3
      t19 = t24 * t7 + t5
      t5 = t15 * t7 + t5
      t15 = t9 * t26 - t29 * t27 * t3 * t2 / 4
      t13 = t9 * t13 - t18 * t27 * t3 * t2 / 4
      t12 = 0.1q1 / t12
      t24 = 0.1q1 / t11
      t4 = 0.1q1 / t4
      t26 = 0.1q1 / t3
      t30 = t18 * t19 + t29 * t5
      t32 = t26 * t31 - t1
      t33 = t26 * t28 - t1
      t34 = t31 * t14
      t35 = t8 * t26
      t36 = t4 * t27
      t37 = t36 * t7
      t38 = t31 + t28
      t39 = t38 * t26
      t40 = t6 * t11
      t36 = t36 * t1 * t3
      t41 = t9 * t4
      t42 = t24 * t8
      t43 = t29 + t18
      t44 = t29 * t31
      t45 = t18 * t28
      t46 = t4 * t43
      t16 = t27 * t7 * (t8 * (t12 * (t26 * (t4 * (t2 * (-t45 - t44) + t6
     # * (-t13 * t18 - t15 * t29)) + t6 * t30 * t24) + t46 * t2 * t1) + 
     #t46 * t7) - t16 * t4 * t26 * (t17 * t18 + t22 * t29) * t20 - t46 *
     # t40)
      t3 = t37 * t12 * (t42 * t7 * (t45 + t44) + (t24 * (t10 * t18 + t14
     # * t29) + t43 * t6) * t3 * t1)
      t17 = -t41 * t11 + t7
      t2 = t21 * (-t7 * t11 * t4 * t26 + t12 * t26 * t17 * t6 + t4 * t12
     # * (t26 * (-t11 * t2 + t19 + t5) - t7) * t1) - t1 * t6 * t12 * t17
     # + t1 * t7 * t12 * (t24 * t7 - t41) * t8
      ret = -3 * t3 + 8 * t2 + 2 * t37 * (-t35 * t30 + t12 * (t18 * (t21
     # + t28) + t29 * (t21 + t31)) * t6 + t12 * (t18 * (t33 * t19 * t8 +
     # t10 * t28) + t29 * (t32 * t5 * t8 + t34)) * t24) + 4 * t12 * (t26
     # * (t40 * (t1 * (t15 + t13) - t31 * t9) * t4 - t1 * t24 * (t19 * t
     #10 + t14 * t5)) + t36 * t42 * t25 * t18 + t35 * t4 * (t1 * (-t15 -
     # t13) + t38 * t9) * t7) + 4 * t12 * (t7 * (t24 * (t10 * t33 + t14 
     #* t32 + t35 * t1 * (t19 + t5)) + t39 * t6) + t1 * t26 * (t4 * (t10
     # * t13 + t14 * t15) + t6 * (-t19 - t5)) + t41 * (t1 * (t14 + t10) 
     #+ t26 * (-t28 * (t40 + t10) - t34)) + t42 * (t36 * t29 - t39) * t2
     #5) + t42 * t23 * t27 ** 2 * t25 * t4 * t12 * (t18 ** 2 + t29 ** 2)
     # / 4 + t16

      hjetmass_qqbgg_box_pmpp_0_s12_0_mhsq_s123_s124 = ret/32q0/(0,1q0)
      return

      end function

      complex*32 function hjetmass_qqbgg_bubble_pmmp_mhsq
     &     (i1,i2,i3,i4,za,zb,mt,p,flip)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          double precision mt
          real*16 p(mxpart,4)
          complex*32 hjetmass_qqbgg_bubble_pmmp_s124
          complex*32 hjetmass_qqbgg_bubble_pmmp_s123
          complex*32 hjetmass_qqbgg_bubble_pmmp_s34
          complex*32 hjetmass_qqbgg_bubble_pmmp_s12
          logical flip



      hjetmass_qqbgg_bubble_pmmp_mhsq = 
     & -hjetmass_qqbgg_bubble_pmmp_s124(i1,i2,i3,i4,za,zb,mt,p,flip)
     & -hjetmass_qqbgg_bubble_pmmp_s123(i1,i2,i3,i4,za,zb,mt,p,flip)
     & -hjetmass_qqbgg_bubble_pmmp_s34(i1,i2,i3,i4,za,zb,mt)
     & -hjetmass_qqbgg_bubble_pmmp_s12(i1,i2,i3,i4,za,zb,mt)
      return

      end function

      complex*32 function hjetmass_qqbgg_bubble_pmmp_s123
     &     (i1,i2,i3,i4,za,zb,mt,p,flip)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt
          real*16 p(mxpart,4)
          complex*32 alpha
          logical flip

      alpha = (za(i1,i2)*zb(i2,i1) + za(i1,i3)*zb(i3,i1) +
     & za(i2,i3)*zb(i3,i2))/(za(i1,i2)*zb(i2,i1) + za(i1,i3)*zb(i3,i1))

      p(5,:) = real(alpha)*p(i1,:)
      p(6,:) = (1q0-real(alpha))*p(i1,:) + p(i2,:) + p(i3,:)
      if (flip .eqv. .true.) then
          call spinoru_qp(6,p,zb,za)
      else
          call spinoru_qp(6,p,za,zb)
      end if

      t1 = za(6, i2)
      t2 = za(i2, i3)
      t3 = zb(i1, 6)
      t4 = zb(i2, 5)
      t5 = za(6, i1)
      t6 = za(5, i2)
      t7 = zb(i2, 6)
      t8 = t5 * t3
      t9 = t6 * t4
      t10 = t1 * t7
      t11 = t10 + t8 - t9
      t12 = 4 * t4
      t13 = t6 * t7
      t11 = t1 * t12 * t13 + t11 ** 2
      t11 = sqrt(t11)
      t14 = t10 + t8 - t9 + t11
      t15 = 0.1q1 / t7
      t16 = (0.1q1 / 0.2q1)
      t17 = t14 * t15
      t18 = t17 * t16
      t19 = -t1 + t18
      t20 = za(6, i3)
      t21 = za(5, i3)
      t22 = 0.1q1 / t6
      t23 = t18 * t21 * t22 - t20
      t24 = zb(i4, i2)
      t25 = za(5, i4)
      t26 = za(i1, i2)
      t27 = za(i3, i4)
      t28 = za(6, i4)
      t29 = 0.1q1 / t25
      t30 = t28 * t29
      t18 = t18 * t22
      t31 = t18 - t30
      t11 = t10 + t8 - t9 - t11
      t32 = t15 * t22
      t33 = t32 * (t14 - t11)
      t34 = 0.1q1 / t14
      t35 = 2
      t36 = t1 * t4
      t37 = t36 * t35
      t38 = t34 * t37 + t18
      t39 = zb(i2, i1)
      t40 = zb(i3, i1)
      t41 = zb(i4, 5)
      t42 = zb(i4, 6)
      t43 = 0.1q1 / t42
      t44 = t41 * t43
      t45 = t44 + t18
      t46 = zb(i3, 5)
      t47 = t20 * t46
      t48 = t28 * t41
      t49 = t48 + t36 + t47
      t50 = zb(i3, 6)
      t51 = t21 * t50
      t52 = t25 * t42
      t53 = t52 + t13 + t51
      t54 = t21 * t46
      t55 = t20 * t50
      t56 = t25 * t41
      t57 = t28 * t42
      t58 = t57 - t54 + t10 + t8 - t9 + t55 - t56
      t58 = t58 ** 2
      t59 = 4 * t49 * t53 + t58
      t59 = sqrt(t59)
      t60 = t57 - t54 + t10 + t8 - t9 + t55 - t56 + t59
      t61 = 0.1q1 / t53
      t62 = t60 * t61
      t63 = t62 * t16
      t64 = t6 * t63 - t1
      t65 = t21 * t63 - t20
      t66 = t63 * t7 + t4
      t67 = zb(i4, i1)
      t58 = 4 * t49 * t53 + t58
      t58 = sqrt(t58)
      t68 = t57 - t54 + t10 + t8 - t9 + t55 - t56 + t58
      t69 = t32 * t11
      t70 = -t69 + t62
      t71 = t17 * t22
      t62 = -t71 + t62
      t68 = 0.1q1 / t68
      t72 = t35 * t49
      t73 = t68 * t72 + t63
      t59 = t57 - t54 + t10 + t8 - t9 + t55 - t56 - t59
      t74 = t61 * (-t60 + t59)
      t75 = t11 * t15
      t76 = t16 * t75 - t1
      t77 = t69 * t16
      t78 = t21 * t77 - t20
      t79 = t77 - t30
      t80 = 0.1q1 / t11
      t37 = t37 * t80 + t77
      t81 = t44 + t77
      t82 = za(i1, i3)
      t83 = t59 * t61
      t84 = t83 * t16
      t85 = t6 * t84 - t1
      t86 = t21 * t84 - t20
      t58 = t57 - t54 + t10 + t8 - t9 + t55 - t56 - t58
      t87 = -t69 + t83
      t83 = -t71 + t83
      t58 = 0.1q1 / t58
      t72 = t72 * t58 + t84
      t88 = t7 * t84 + t4
      t89 = t21 * t28
      t90 = t29 * t89 - t20
      t91 = -t28 * t90
      t92 = zb(i3, i2)
      t93 = t26 * t39
      t94 = t82 * t40
      t95 = t2 * t92
      t96 = t95 + t93 + t94
      t97 = t44 + t30
      t98 = t44 * t6 + t1
      t77 = t42 * t77 + t41
      t99 = zb(i4, i3)
      t100 = za(i2, i4) * t24
      t101 = t27 * t99
      t102 = t67 * za(i1, i4)
      t103 = -t102 - t94 - t95 - t55 + t57 - t100 - t101
      t104 = t4 * t21
      t105 = t7 * t20
      t106 = t105 + t104
      t107 = t20 * t4
      t108 = t107 * t42
      t109 = t8 + t93
      t110 = t36 * t21
      t111 = t55 + t8
      t112 = t1 * t28
      t113 = t36 * t20
      t114 = t15 ** 2
      t115 = t15 * t114
      t116 = t11 ** 2
      t117 = t116 ** 2
      t118 = t11 * t116
      t119 = t1 ** 2
      t120 = t1 * t119
      t121 = t22 ** 2
      t122 = t119 * t7
      t123 = t110 * t15
      t124 = t57 * t8
      t125 = t93 * t111
      t126 = t1 * (t112 * (t106 * t41 + t108) + t113 * (t102 + t94 + t95
     # + t100 + t101) + t47 * (t109 * t20 - t110))
      t127 = t115 * t21
      t128 = t11 * t22
      t129 = t128 * t1
      t130 = t89 * t41
      t131 = t20 * (t102 + t94 + t95 - t54 + t57 + t100 + t101) + t130
      t132 = t25 * t20
      t133 = -t89 - t132
      t134 = t57 * t20
      t135 = t133 * t41 + t134
      t136 = t4 ** 2
      t137 = t20 ** 2
      t138 = t6 * t46
      t139 = t119 * t50
      t140 = t137 * t46
      t141 = t21 * t120
      t142 = t8 * t69
      t143 = t119 * t121
      t144 = t143 * t21
      t145 = t21 * t116 * t114 * t22
      t146 = t20 * t1
      t147 = t146 + t145
      t148 = t132 * t41
      t149 = -t148 + t110
      t150 = t3 ** 2
      t151 = t1 * t21
      t152 = t6 * t28
      t153 = t100 * t1
      t154 = t101 * t1
      t155 = t22 * t4
      t156 = t155 * t127
      t157 = t20 * t22
      t158 = t157 * t114
      t159 = t20 * t41
      t160 = t159 * t112
      t161 = t150 * t5 ** 2
      t145 = t11 * (t15 * (t116 * (t15 * (t149 * t22 + t151 * (t102 - t5
     #7 + t95 + t55 + t94) * t121) - t108 * t25 * t22 * t114) + t152 * t
     #107 * t41) + t153 * t22 * t147 + t154 * t22 * t147) + t126 + t129 
     #* (t20 * (t109 * t15 * t55 + t122) + t123 * t103) + t127 * (-t124 
     #+ t125) * t121 * t118 + t128 * t119 * (t20 * (t102 + t57 - t54 + t
     #95 + t55 + t94) + t130) + t137 * (t138 * t75 + t139) * t4 + t141 *
     # t136 - t144 * t118 * t15 + t93 * t75 * (t1 * t135 * t22 - t140) +
     # t140 * t122 + t142 * t1 * t131 + t8 * (t118 * (t158 + t156) + t16
     #0) + t93 * (t118 * (t158 - t156) + t160 - t159 * t75 * t28) + t161
     # * t69 * (t146 - t145)
      t147 = -t89 - t132
      t156 = t20 * t42
      t162 = t21 * t41
      t163 = -t162 + t156
      t164 = t48 + t47
      t165 = t21 ** 2
      t166 = t165 * t46
      t167 = t147 * t41
      t168 = t141 * t121
      t169 = t151 * t136
      t170 = t20 * t15
      t171 = t8 * t22
      t172 = t171 * t114
      t173 = t93 * t22
      t174 = t173 * t114
      t175 = t107 * t15
      t176 = -t175 + t144
      t177 = -t57 + t55
      t178 = t95 + t94
      t179 = t89 * t1
      t180 = t179 * t42 * t22
      t181 = t137 * t50
      t182 = t161 * t151
      t183 = t182 * t121
      t184 = t183 * t15
      t185 = t94 * t176
      t186 = t93 * t32
      t187 = t15 * t116
      t188 = t102 + t94 + t95 + t55 - t57 + t100 + t101
      t189 = t56 + t9
      t190 = t151 * t22
      t191 = t190 + t20
      t192 = t1 * t50
      t193 = t192 * t22
      t194 = t21 * t22
      t195 = t132 * t1
      t196 = t195 * t42
      t197 = t175 * t189
      t198 = t8 * t15
      t199 = t93 * t15
      t200 = t8 - t93
      t201 = t132 * t42
      t202 = t200 * t21
      t203 = t202 - t201
      t204 = t203 * t22 * t1
      t205 = t32 * t116
      t206 = t205 * t1 * (t20 * (-t102 - t94 + t54 - t57 - t100 - t101) 
     #+ t204) + t187 * (t199 * (t20 * (t194 * (-t193 + t46) + t4) - t171
     # * t191) - t197 + (t20 * (-t55 - t95) - t130 + t104 * t188 * t15) 
     #* t22 * t1 + t198 * (t22 * (-t196 * t22 - t130) + t107))
      t207 = t162 * t25
      t208 = (t89 + t132) * t42
      t209 = t208 - t207
      t210 = t52 + t51
      t211 = t93 * t209
      t212 = t114 * t22
      t213 = t212 * t118
      t214 = -t102 - t94 - t95 - t55 + t57 - t100 - t101
      t215 = t166 * t32
      t216 = t213 * (t104 * t15 * t214 + t200 * t215 + t134) + t213 * (t
     #15 * (t22 * (t8 * (t21 * (t102 + t94 + t95 - t55 + t56 + t100 + t1
     #01) + t201) + t211) + t151 * t93 * t210 * t121) + t20 * (t102 + t9
     #4 + t95 - t54 + t55 + t100 + t101) + t130)
      t217 = t104 * t15
      t218 = t217 + t20
      t219 = t4 * t42
      t220 = t15 * t219 + t41
      t221 = -t217 + t20
      t222 = t198 * t218
      t223 = t132 * t15 * t220
      t224 = t48 + t36
      t225 = t128 + t4
      t226 = t9 + t11
      t227 = t1 * t25
      t228 = t8 * t20
      t229 = t4 * t50
      t230 = t119 * t22
      t231 = t8 * t1
      t232 = t93 * t1
      t233 = t102 + t94 + t95 - t54 + t57 + t100 + t101
      t234 = t75 * t1
      t149 = t234 * (t121 * t187 * t203 + t107 * t233) + t20 * (t219 * t
     #69 * t119 * t25 + t4 * t118 * t114 - t120 * t4 * t7 + t1 * (t225 *
     # t227 + t226 * t28) * t41) + t1 * (t226 * t46 + t229 * t75) * t137
     # + t231 * (t128 * (-t149 * t15 - t146) - t113) + t110 * t11 * (t15
     # * t224 - t230) + t232 * (t128 * (t15 * (t21 * (t47 + t36) - t228)
     # - t146) - t113)
      t226 = t8 * t210
      t235 = t93 * t210
      t236 = t15 * (t32 * (-t235 + t226) - t94 - t95 - t55 + t57 - t100 
     #- t101 - t102) + t1
      t237 = t121 * t114
      t238 = t237 * t117
      t239 = -t56 + t93 + t8
      t240 = t36 * t6
      t241 = -t8 + t93
      t242 = t241 * t21
      t243 = t242 + t201
      t244 = (0.5q1 / 0.8q1)
      t245 = -(0.5q1 / 0.16q2)
      t246 = (0.3q1 / 0.2q1)
      t247 = (0.3q1 / 0.4q1)
      t248 = -(0.15q2 / 0.2q1) * t4
      t249 = t136 * t6
      t250 = 5 * t249
      t251 = 8 * t36
      t252 = t20 * t11
      t253 = t54 - t8
      t254 = t55 * t1
      t255 = t192 * t137
      t256 = t120 * t20 * t121
      t257 = t93 * t21
      t233 = t233 * t20
      t258 = (t167 + t134) * t22
      t259 = t198 * t121 * t1
      t260 = t175 - t144
      t261 = t69 * t25
      t262 = -t261 + t28
      t263 = -t75 + t1
      t264 = t151 * t75 * t121
      t261 = (t261 - t28) * t42
      t265 = t93 * t69
      t158 = t187 * (t170 * t161 * t1 * t121 + t143 * (t20 * (t262 * t42
     # + t55) + t130) + t142 * (t175 - t264 + t144) + t175 * t41 * t262 
     #+ t146 * t102 * t121 * t263 + t146 * t94 * t121 * t263 + t146 * t9
     #5 * t121 * t263 + t234 * t121 * (t20 * (-t55 + t261) - t130) + t26
     #5 * t260) + t187 * (t20 * (t175 * t46 + t143 * (t101 - t54 + t100)
     # + t234 * (-t101 + t54 - t100) * t121) + t186 * (t1 * (t75 * t21 *
     # (-t55 + t11) * t121 + t258) - t140) + t259 * (t69 * (-t201 - t257
     #) + t233 + t130)) + t116 * (t4 * (t144 * t214 * t114 + t128 * (t18
     #8 * t190 - t252) * t115) + t256 + t158 * t93 * (t22 * (t253 * t75 
     #+ t254) - t48) - t252 * t136 * t115 + t237 * t8 * (-t130 * t75 + t
     #255))
      t234 = t151 * t15
      t266 = t8 * t32
      t267 = t121 * t115
      t268 = t267 * t117
      t269 = t217 - t20
      t270 = -t190 - t20
      t271 = (t193 + t46) * t21
      t272 = t89 * t22
      t273 = t151 * t46
      t274 = t140 * t15
      t264 = t175 + t264 - t144
      t275 = t75 - t1
      t276 = -t8 + t11
      t277 = t1 * t15
      t278 = t132 * t69
      t279 = t237 * t161
      t280 = t162 + t156
      t281 = t280 * t28
      t264 = t213 * (t104 * t114 * t128 * t8 + t143 * t42 * t89 + t15 * 
     #t281 * t4 - t181 * t186 + t264 * t94 + t264 * t95 + t184) + t213 *
     # (t20 * (t175 * t50 + t48) + t21 * (t254 * t121 * t275 + t15 * (t1
     #29 - t47) * t4 + t277 * t136 - t112 * t75 * t42 * t121) + t102 * (
     #t121 * t234 * t276 - t144 + t175) + t101 * t264 + t100 * t264 - t2
     #78 * t220 + t186 * (t11 * t221 - t180) + t142 * t20 - t279 * t21 *
     # t11) + t32 * t118 * (-t144 * t75 + t168 + t274 + t174 * (t22 * (t
     #21 * (t227 * t41 + t55 * t75 + t273) - t196) + t28 * (t162 - t156)
     # + t148) + t172 * (t20 * (-t102 + t56 - t271) + t257 * t69 + t272 
     #* t263 * t42 + t100 * t270 + t94 * t270 + t101 * t270 + t95 * t270
     #))
      t282 = t8 + t93
      t283 = t75 * t4
      t284 = t22 * t227 + t28
      t285 = t128 + t4
      t286 = t32 * t1
      t287 = t41 * t284
      t288 = t170 * t93
      t123 = t129 * (t137 * (t15 * t229 * t263 + t263 * t46) + t20 * (t1
     #1 * (-t101 * t114 * t4 - t15 * t287) + t116 * (t212 * t56 - t286) 
     #+ t112 * t41) + t170 * t8 * (t128 * (t15 * (-t56 + t93) + t1) + t4
     #8 - t205) + t123 * (t1 * t285 - t285 * t75) + t288 * (t128 * (-t15
     # * (t54 + t11) + t1) + t48)) + t69 * t1 * (t155 * t243 * t114 * t1
     #16 + t36 * (t130 + t233) + t140 * t282 + t283 * (t20 * (-t102 - t9
     #4 - t95 + t54 - t57 - t100) - t130 + t204))
      t204 = t119 * t4
      t205 = t127 * t22 * t121
      t233 = -(0.5q1 / 0.4q1)
      t285 = (0.1q1 / 0.4q1)
      t289 = (0.1q1 / 0.8q1)
      t290 = (0.1q1 / 0.16q2)
      t291 = (0.1q1 / 0.32q2)
      t123 = -t158 * t285 - t16 * t123 + t35 * t107 * t75 * t119 * t225 
     #+ t290 * (t268 * (t20 * (-t101 + t261 - t100) + t21 * ((-t57 + t94
     #) * t15 * t4 - t48) + t102 * t269 + t95 * t269) + t268 * (t20 * (-
     #t55 + t11 - t94 + t54) + t93 * (-t234 * t210 * t121 + t32 * (t21 *
     # (-t57 + t11 + t54 + t56) - t201)) + t217 * (t101 + t55 + t100) + 
     #t266 * (t21 * (-t102 - t11 - t94 - t95 - t54 + t55 - t56 - t100 - 
     #t101) - t201))) - t289 * t264 - t291 * t205 * t11 * t117 * t236 + 
     #t233 * t113 * t116 * t114 * t225 + t113 * (t142 * t275 - t204 + t1
     #28 * (t56 * (-t11 * t114 + t277) - t119) + t48 * t275 + t47 * t275
     # + t265 * t275) + t146 * t238 * t245
      t18 = t18 * t42 + t41
      t142 = t152 * t29
      t158 = t142 - t1
      t84 = t42 * t84 + t41
      t63 = t42 * t63 + t41
      t234 = 0.1q1 / t53
      t261 = t228 - t110
      t264 = t6 * t20
      t265 = -t264 - t151
      t268 = t156 * t112
      t275 = t234 ** 2
      t292 = t275 ** 2
      t293 = t234 * t275
      t294 = t59 ** 2
      t295 = t294 ** 2
      t296 = t59 * t294
      t297 = t93 * (t41 * (t265 * t28 - t195) + t268)
      t298 = t6 * t21
      t299 = t105 - t104
      t300 = t7 ** 2
      t301 = t6 ** 2
      t302 = t119 * t300
      t303 = t107 * t301
      t304 = t122 * t21
      t305 = t94 * t119
      t306 = t305 * t299
      t307 = t95 * t119
      t308 = t93 * t137
      t309 = t48 * t7
      t310 = t7 * t41
      t311 = t310 + t219
      t312 = -t10 + t9
      t313 = t4 * t301
      t314 = t313 * t46
      t315 = t122 * t50
      t316 = t101 * t119 * t299
      t317 = t102 * t119
      t318 = t100 * t119
      t319 = t8 * t6
      t320 = t93 * t6
      t321 = t151 * t7
      t322 = t161 * t21
      t323 = t56 - t93 - t8
      t324 = t234 * t59
      t325 = t113 * t6
      t326 = t48 + t47
      t327 = t7 * t46
      t328 = t327 + t229
      t329 = t8 - t93
      t330 = t227 * t311
      t331 = t6 * t137
      t332 = t151 * t329
      t333 = t8 * t41
      t334 = t324 * t146
      t335 = t334 * (t25 * (-t13 * t275 * t294 * t42 + t333) + t9 * (t10
     #2 + t94 + t95 - t54 + t100 + t101) + t93 * t253) + t113 * (t1 * (t
     #56 - t93 - t8) + t326 * t6 - t122) + t324 * t1 * (t20 * (t152 * t3
     #11 + t330) + t104 * (t224 * t6 - t122) - t231 * t106 - t232 * t299
     # + t331 * t328) + t13 * (t332 + t303) * t293 * t296
      t336 = t310 + t219
      t337 = t105 - t104
      t338 = t105 + t104
      t339 = t93 * t337
      t340 = t8 * t338
      t341 = t10 * t338 - t132 * t336 + t339 + t340
      t342 = t294 * t275
      t343 = t342 * t6
      t344 = t227 * t42
      t345 = t264 + t151
      t346 = t298 * t47
      t313 = (-t313 + t344) * t20
      t347 = (t48 - t47) * t21
      t348 = t154 * t6
      t349 = t348 * t337
      t350 = t102 * t1
      t351 = t350 * t6
      t352 = t351 * t337
      t353 = t153 * t6
      t354 = t353 * t337
      t355 = -t9 + t93
      t356 = t338 * t42
      t357 = t356 * t1
      t358 = t94 * t1
      t359 = t358 * t6
      t360 = t95 * t1 * t6
      t361 = t151 * t55
      t362 = t342 * (t361 * t355 + t152 * (t162 * t8 + t357) + t359 * t3
     #37 + t360 * t337) + t342 * (t13 * (t347 + t181) * t1 + t349 + t352
     # + t354 + t93 * (t304 - t303 - t346) + t8 * (t345 * t93 - t304 + t
     #313) + t201 * t122 + t303 * t189)
      t363 = -t105 - t104
      t364 = t363 * t41 - t108
      t365 = t54 * t107
      t366 = t304 - t303
      t367 = -t264 - t151
      t368 = t166 * t1
      t369 = t181 * t6
      t141 = t141 * t300
      t370 = t95 * t366
      t371 = t101 * t366
      t372 = t102 * t366
      t373 = t100 * t366
      t374 = t304 * t177
      t375 = (t28 * t364 + t365) * t301
      t376 = -t327 - t229
      t208 = t342 * (t301 * (t137 * t376 - t169) - t182 + t93 * (t1 * t2
     #08 + t152 * t163) + t8 * (t151 * t188 + t264 * (t102 + t94 + t95 -
     # t56 + t100 + t101))) + t342 * (t94 * t366 - t141 + t370 + t371 + 
     #t372 + t373 + t374 + t375 + t93 * (t367 * t56 - t368 + t369) + t34
     #6 * t8)
      t366 = t52 + t13
      t377 = t366 * t20
      t378 = t292 * t295
      t379 = t13 * t20
      t380 = t21 * (-t9 + t8) + t379
      t381 = t57 - t55
      t382 = -t138 + t192
      t211 = t6 * t296 * t293 * (t8 * (t21 * (t102 + t94 + t95 - t55 + t
     #56 + t100) + t201) + t211) + t296 * (t6 * (t293 * (t100 * t6 * t33
     #7 + t102 * t6 * t337 + t104 * t6 * t381 + t94 * t6 * t337 + t13 * 
     #(t20 * (t57 - t54 + t55) + t130) + t166 * t8) + t95 * t6 * t293 * 
     #t337 + t101 * t293 * t380) + t257 * t293 * (t21 * t382 + t344))
      t383 = t1 * t300
      t384 = t6 * (t214 * t7 + t383) + t226 - t235
      t385 = t4 * (-5 * t9 + 8 * t10)
      t334 = t334 * t6 * (t248 * t324 * t13 + t385 - (0.5q1 / 0.2q1) * t
     #343 * t300)
      t386 = t48 + t36 + t47
      t387 = t138 - t192
      t388 = -t105 + t104
      t389 = t363 * t42
      t390 = t257 * t46
      t391 = t136 * t301
      t392 = t130 * t6
      t393 = t146 * t7
      t327 = t327 + t229
      t394 = t55 + t8
      t395 = t6 * t1
      t396 = t146 * t161
      t397 = t137 * (t315 + t314)
      t255 = t342 * (t231 * (t20 * (-t57 - t94 - t95 + t54) - t130) - t3
     #96 - t397 + t93 * (t41 * (t28 * t345 + t195) - t268) + t342 * t107
     # * t6 * t301 * t7 + t324 * (t93 * (t151 * t394 + t228 * t6) + t395
     # * (-t21 * t20 * t327 + t94 * t337 + t95 * t337))) + t342 * (t1 * 
     #(-t228 * (t102 + t100) + t101 * (t20 * (-t10 - t8) + t110) - t10 *
     # (t130 + t393)) + t342 * t13 * t1 * t203 + t324 * (t6 * (t20 * (t3
     #91 - t390) + t112 * (t162 * t7 + t356)) + t349 + t352 + t354 + t8 
     #* (t392 + t196))) + t342 * (t28 * (t119 * t389 - t303 * t41) + t30
     #8 * t387 + t305 * t388 + t307 * t388 + t317 * t388 + t318 * t388 +
     # t304 * t47 + t324 * (t122 * t243 + t13 * t255 + t303 * t323) + t2
     #54 * (t110 - t228))
      t108 = t338 * t41 + t108
      t349 = -t304 + t303
      t219 = -t310 - t219
      t310 = t110 * t7
      t352 = t93 * (-t369 + t368)
      t354 = t95 * t349
      t356 = t101 * t349
      t398 = t94 * t349
      t399 = t102 * t349
      t349 = t100 * t349
      t365 = (t108 * t28 - t365) * t301
      t400 = t304 * t381
      t401 = t324 * t301
      t402 = t293 * t296
      t403 = t89 * t6
      t404 = t28 * t367 - t195
      t405 = -t102 - t100 - t101
      t406 = t93 * (t41 * (t25 * t345 + t403) + t42 * t404)
      t407 = t331 * t327
      t408 = t169 * t6
      t409 = t324 * t6
      t410 = t47 * t21
      t411 = (-t162 - t156) * t28
      t412 = t101 * t388
      t413 = t94 * t388
      t414 = t102 * t388
      t415 = t100 * t388
      t416 = t377 + t257
      t387 = t93 * t165 * t387
      t417 = t298 * t8
      t418 = t301 * (t7 * (t411 + t410) + t412 + t413 + t414 + t415 + t9
     #5 * t388 + t104 * t177)
      t419 = t93 * (t42 * (t25 * t367 - t403) + t207 * t6)
      t420 = t181 * t301 * t7
      t421 = t47 * t7
      t422 = t232 * t388
      t423 = t231 * t338
      t364 = t61 * (t112 * t364 + t36 * (-t20 * (t102 + t94 + t95 + t100
     # + t101) - t110) + t47 * (t20 * (-t8 - t93) + t110))
      t424 = t110 * t13
      t425 = t59 * t1
      t426 = t104 * t122
      t427 = t228 * (t56 - t93)
      t282 = t61 * t20 * (t146 * t376 - t48 * t282)
      t376 = t6 * (t20 * (t219 * t25 + t383) + t339 + t340) * t293
      t294 = t425 * (t376 * t294 + t282 + (t20 * (t6 * (t4 * (t100 + t94
     # + t95 + t57) + t309) + t54 * t355) + t408 - t426 + t427 + t9 * t1
     #30) * t275 * t59) + t425 * (t424 * t294 * t293 + t364 + (t20 * (t6
     # * (t4 * (t102 + t55 + t101) + t421) + t227 * t336) + t422 - t423)
     # * t275 * t59)
      t355 = t324 * t7 + t4
      t226 = t6 * (t188 * t7 - t383) - t226 + t235
      t235 = t298 * t234 * t292
      t383 = t409 * t107 * t119 * t355
      t425 = t378 * t146 * t245 * t301 * t300
      t188 = t16 * t294 + t233 * t342 * t113 * t301 * t355 + t285 * t255
     # - t289 * (t402 * (t6 * (-t302 * t324 * t21 + t407 + t408) + t8 * 
     #(-t151 * (t100 + t95) + t264 * (t56 - t94 - t95 - t54)) + t322 * (
     #-t409 + t1)) + t402 * (t324 * t298 * (t10 * t188 + t93 * t394 - t1
     #24) + t141 + t406 + t8 * (t264 * t405 + t151 * (-t102 - t94 - t55 
     #+ t57 - t101))) + t402 * (t401 * (t132 * t219 + t310 + t339 + t340
     #) + t352 + t354 + t356 + t398 + t399 + t349 + t365 + t400)) + t290
     # * (t378 * (t319 * (t21 * (-t102 - t94 - t95 + t55 - t56) - t201) 
     #+ t419 - t420) + t378 * (t401 * t7 * t416 + t387 + t417 * (-t324 *
     # t13 - t100 - t101 - t54) + t418)) + t291 * t235 * t59 * t295 * t2
     #26 + t35 * t383 - t113 * (t1 * t386 + t324 * (t1 * t239 + t6 * (-t
     #48 - t47) + t122) + t343 * t323) + t425
      t219 = t60 ** 2
      t255 = t219 ** 2
      t294 = t60 * t219
      t295 = t48 * t301
      t139 = t139 * t21
      t339 = t21 * t119
      t340 = t159 * t28
      t355 = t9 - t93
      t120 = t120 * t300
      t383 = t47 * t6
      t401 = t231 * t20
      t402 = t1 * (t340 * t109 + t339 * t136)
      t229 = t1 * (t201 * t13 * t294 * t293 + t107 * (t1 * (-t56 + t93 +
     # t8 + t10) - t383) + (t6 * (t4 * (t20 * (-t102 - t94 - t95 + t54 -
     # t57 - t101) - t130) - t105 * t164) + t93 * (t21 * (-t47 - t36) + 
     #t393) + t8 * (t20 * (-t56 + t93 + t10) + t110)) * t234 * t60) - t1
     #52 * t113 * t41 + t13 * (-t332 - t303) * t293 * t294 + t1 * (t20 *
     # (-t9 * t100 - t330) - t331 * t229 - t110 * t312) * t234 * t60
      t330 = t294 * (t293 * (-t420 + t419) + t319 * t293 * (t21 * (-t102
     # - t94 + t55 - t100) - t201)) + t294 * (t387 * t293 + t418 * t293 
     #+ t417 * t293 * (-t101 - t95 - t54 - t56))
      t292 = t292 * t255
      t331 = t8 * t21
      t332 = t275 * t219
      t182 = t332 * (t169 * t301 + t141 + t182 + t406 + t8 * (t264 * (-t
     #102 - t94 - t95 + t56 - t100 - t101) + t151 * (-t102 - t95 + t57 -
     # t100)) + t301 * t137 * t327) + t332 * (t331 * (t1 * (-t101 - t94 
     #- t55) - t383) + t352 + t354 + t356 + t398 + t399 + t349 + t365 + 
     #t400)
      t327 = t8 * t367
      t349 = t151 * t50
      t352 = -t56 * t6 + t349
      t354 = t107 * t6
      t356 = t359 * t388
      t346 = t93 * (-t361 + t303 + t346 + t327 - t304)
      t359 = t391 * t20
      t360 = t360 * t388
      t348 = t332 * (t6 * (t389 * t112 - t359) + t360 + t351 * t388 + t2
     #31 * (t321 - t201)) + t332 * (t354 * t352 + t346 + t356 + t10 * (t
     #20 * (t138 * t21 - t55 * t6 - t344) - t392) + t353 * t388 + t348 *
     # t388 + t319 * (t354 - t130))
      t351 = t234 * t60
      t353 = t325 * (t351 * t239 + t36)
      t361 = t351 * t146 * t6 * (t351 * t248 * t13 - (0.5q1 / 0.2q1) * t
     #332 * t6 * t300 + t385)
      t365 = t10 + t8
      t385 = -t55 + t54
      t387 = t102 + t55 + t100 + t101
      t389 = t351 * t6
      t391 = t1 * t60
      t108 = t391 * (t424 * t219 * t293 + t364 + (t20 * (t6 * (t387 * t4
     # + t421) + t52 * t36) + t422 - t423) * t275 * t60) + t391 * (t376 
     #* t219 + t282 + (t1 * (t148 * t7 + t249 * t21) + t20 * (t9 * (-t54
     # + t94 + t95) + t390) - t426 + t427 + t152 * t108) * t275 * t60)
      t219 = t113 * (t1 * (-t48 - t36 - t47) + t351 * (t1 * t323 + t164 
     #* t6 - t122) + t332 * t239 * t6)
      t282 = t389 - t1
      t364 = t351 * t7
      t376 = t293 * t294
      t390 = t8 * (t102 + t100 + t101)
      t391 = t351 * t13
      t141 = t376 * (t264 * (t389 * t336 * t25 + t390) - t141 + t151 * (
     #t391 * t405 + t390)) + t376 * (t389 * (-t319 * t338 + t321 * (t10 
     #- t9) + t320 * t388) + t370 + t371 + t372 + t373 + t374 + t375 + t
     #93 * (t369 - t368) + t94 * (t151 * t365 - t303)) + t376 * (t6 * (t
     #151 * (t364 * (t57 - t94 - t95 - t55) - t249) - t407) + t93 * (t1 
     #* t209 + t6 * (t20 * (-t51 * t351 - t56 + t57) - t130)) + t322 * t
     #282 + t8 * (t20 * (t6 * (-t56 + t94 + t54) + t349) + t95 * t345 + 
     #t89 * t282 * t42 - t389 * t257))
      t249 = t364 + t4
      t282 = t292 * t6 * (t6 * (t364 * t416 - t105 * (t57 + t95) - t89 *
     # t336) - t331 * (t391 + t100 + t101 + t102)) + t292 * (t301 * (t10
     #5 * t385 + t412 + t413 + t414 + t415 + t104 * (t55 + t95)) + t93 *
     # (t6 * (t21 * (-t57 + t54 + t56) - t201) - t151 * t210) + t319 * (
     #t21 * (-t56 - t94 - t95 - t54 + t55) - t201))
      t108 = t16 * t108 + t233 * t332 * t113 * t301 * t249 - t285 * (t33
     #2 * (t20 * (t4 * (-t139 + t295) + t120 + t122 * (t101 - t54)) + t3
     #05 * t337 + t307 * t337 + t308 * t382 + t389 * (t153 * t388 + t154
     # * t388 + t350 * t388 - t359) + t401 * t387) + t332 * (t1 * (t1 * 
     #(t281 * t7 + t104 * (-t101 + t57)) + t153 * t337 + t350 * t337 + t
     #228 * t178) + t351 * (t20 * (t9 * t352 + t10 * (t385 * t6 - t344))
     # + t346 + t8 * (t304 + t303))) + t332 * (t332 * t13 * (t151 * t241
     # + t313) + t396 + t397 + t93 * (t404 * t41 + t268) + t351 * (-t196
     # * t8 + t356 + t360 + t152 * (-t162 * t365 - t357)) + t231 * (t281
     # - t410))) + t289 * t141 + t290 * t282 + t291 * t235 * t60 * t255 
     #* t226 + t35 * t351 * t354 * t119 * t249 + t219 + t292 * t146 * t2
     #45 * t301 * t300
      t141 = -t102 - t93 - t94 - t95 + t8 + t10 - t55 + t57 - t100 - t10
     #1
      t154 = (t102 + t93 + t94 + t95 + t8 - t54 + t55 + t100 + t101) * t
     #20
      t156 = (t162 + t156) * t28
      t196 = t106 * t1
      t219 = t377 * t275
      t226 = -t16 * t59 * (t61 * (t1 * t338 + (t102 + t93 + t94 + t95 + 
     #t8 - t54 + t55 + t100 + t101) * t20 + t281) - t219 * t59) - t285 *
     # t342 * t21 * (-t102 - t93 - t94 - t95 + t8 + t10 - t55 + t57 - t1
     #00 - t101) - t20 * (t324 * (-t56 - t9) + t36 + t47 + t48)
      t219 = t16 * t60 * (t61 * (t1 * t363 + t20 * (-t102 - t93 - t94 - 
     #t95 - t8 + t54 - t55 - t100 - t101) + t411) + t60 * t219) + t285 *
     # t332 * t21 * (t102 + t93 + t94 + t95 - t8 - t10 + t55 - t57 + t10
     #0 + t101) + t20 * (t351 * t189 - t36 - t47 - t48)
      t235 = t29 ** 2
      t241 = t28 ** 2
      t249 = t28 * t241
      t255 = t107 * t119
      t268 = t241 * t21 * t7 * t235
      t281 = -t35 * t395 * t29 * t241 * (t30 * t337 + t107 - t268) + t28
     # * (t6 * t380 * t29 * t235 * t249 + t255 + (-t304 + t303 + t327) *
     # t235 * t241 - t301 * t21 * t241 ** 2 * t7 * t235 ** 2 + t112 * (t
     #20 * t365 - t110) * t29)
      t282 = t29 * t241
      t300 = t298 * t241 * t235
      t313 = t28 * (t30 * t345 - t146 - t300)
      t327 = t44 * t7 - t4
      t331 = t30 * t7 + t4
      t165 = t192 * t165
      t336 = t264 * t293 * t366 + t293 * (t6 * (t102 + t93 + t94 + t95 -
     # t8 - t57 + t100 + t101) + t344) * t21 + t165 * t293
      t61 = (t102 + t93 + t94 + t95 + t8 + t10 + t55 + t57 + t100 + t101
     #) * t61
      t337 = (-t52 - t13 - t51) * t275
      t162 = t100 * t367 + t101 * t367 + t102 * t367 + t264 * t189 + t93
     # * t367 + t94 * t367 + t95 * t367 + t368 + t162 * (-t152 + t227)
      t189 = t289 * t294 * t336 + t16 * t146 * t60 * (t337 * t60 + t61) 
     #+ t285 * t332 * t162 + t146 * (t351 * (-t56 - t9 - t54) + t36 + t4
     #7 + t48)
      t61 = t289 * t296 * t336 + t16 * t146 * t59 * (t337 * t59 + t61) +
     # t285 * t342 * t162 - t146 * (t324 * (t56 + t9 + t54) - t36 - t47 
     #- t48)
      t162 = t32 * t25 * t200 + t28
      t200 = t194 * t114
      t336 = t122 * t28
      t337 = t20 * t28
      t133 = t133 * t22 * t1
      t338 = t112 * t20
      t345 = t89 * t119
      t346 = t8 * t112
      t349 = t198 + t1
      t350 = t175 * t25 + t272 * t349
      t352 = t89 * t4
      t356 = t186 * (t191 * t25 + t89)
      t357 = t15 * (t171 * t132 + t352) + t337 + t356
      t359 = t186 * (t147 * t22 * t1 - t337)
      t360 = t89 * t121
      t363 = t28 * t218
      t364 = t132 * t32
      t127 = t127 * t121
      t369 = t212 * t25
      t129 = -t16 * t129 * (t28 * (t20 * (t15 * (t8 - t11 + t93) + t1) +
     # t217 * t263) + t278 * (t15 * t276 - t1)) - t285 * t187 * (t175 * 
     #t262 + t359 + t360 * (t1 * t263 + t198 * t263)) - t289 * t213 * (t
     #364 * (t8 - t11) + t356 + t363) - t290 * t127 * t117 * t162 + t113
     # * (-t369 * t116 + t75 * t284 - t112)
      t262 = t14 * t22
      t263 = t14 ** 2
      t276 = t263 ** 2
      t278 = t14 * t263
      t370 = t15 * t263
      t371 = t262 + t4
      t372 = -t17 + t1
      t373 = t217 * t372
      t374 = t262 * t1
      t375 = t71 * t25
      t376 = -t375 + t28
      t380 = t212 * t278
      t382 = t369 * t263
      t284 = -t16 * t374 * (t28 * (t20 * (t15 * (t8 - t14 + t93) + t1) +
     # t373) + t132 * t71 * (t15 * (-t8 + t14) - t1)) - t285 * t370 * (t
     #175 * t376 + t359 + t360 * (t1 * t372 + t198 * t372)) - t289 * t38
     #0 * (t364 * (t8 - t14) + t356 + t363) - t290 * t127 * t276 * t162 
     #+ t113 * (t17 * t284 - t112 - t382)
      t356 = t49 * t20
      t359 = t128 * t16 + t4
      t360 = t262 * t16 + t4
      t363 = t237 * t21
      t364 = t169 * t15
      t385 = t100 + t94 + t95
      t387 = t102 + t101
      t388 = t266 * t21
      t171 = t22 * t115 * t278 * (t171 * (t21 * (t101 + t56 + t100) + t2
     #01) + t93 * (t52 * t151 * t121 + t209 * t22) - t104 * (t55 + t94))
     # + t380 * (t20 * (t102 + t94 + t95 - t54 + t55 + t57 + t100 + t101
     #) + t130 + t93 * (t165 * t15 * t121 - t215) + t388 * (t102 + t94 +
     # t95 + t54 - t55) + t217 * (-t102 - t95 + t57 - t100 - t101))
      t209 = t363 * t278
      t215 = t15 * t4
      t389 = t36 * t42
      t194 = t20 * (t14 * (t1 * (t215 * t385 + t47) + t4 * t114 * t263) 
     #+ t152 * t36 * t41 + t232 * t71 * t253 + t227 * t71 * (t333 + t389
     #)) + t1 * (t1 * (t17 * t136 * t21 + t148 * (t262 + t4)) + t20 * (-
     #t52 * t237 * t278 + t17 * (t102 - t54 + t55 + t101) * t4 + t138 * 
     #t107) - t204 * (t194 * t14 + t105) + t14 * (t280 * t15 * t4 + t159
     #) * t28 + t93 * (t374 * t269 - t113 - t209) + t8 * (-t374 * t218 -
     # t113 + t209))
      t209 = -t217 + t20
      t280 = t71 * t1
      t333 = t339 * t155
      t135 = t17 * (-t333 * (t102 + t101) + (t8 * t131 + t93 * t135) * t
     #22 * t1 + t340 * t355)
      t111 = t278 * (-t267 * t8 * t89 * t42 + t237 * t151 * (t102 + t55 
     #+ t94)) + t402 + t1 * (t393 * t326 + t36 * (t20 * (t102 + t94 + t9
     #5 - t54 + t55 + t57 + t100 + t101) + t130) + t140 * t93) + t278 * 
     #(t114 * (t22 * (t228 + t110) + t151 * (-t57 + t95 + t101 + t100) *
     # t121) + t115 * (-t322 * t121 + t155 * (t21 * t329 - t201))) + t14
     # * (t15 * (t138 * t137 * t4 + t339 * (-t107 * t50 * t22 - t121 * t
     #263)) + t156 * t230 + t157 * t119 * (t102 - t54 + t10 + t55 + t101
     #) + t307 * t22 * t209 + t318 * t22 * t209 + t305 * t22 * t209 - t2
     #12 * t148 * t263) + t231 * t137 * (t71 * t50 + t46) + t93 * t17 * 
     #(t137 * (-t46 + t193) + t263 * (t363 * t111 + t32 * t20)) + t280 *
     # (t161 * t20 + t389 * t89) + t135
      t135 = t286 * t263
      t193 = t135 * (t20 * (t52 * t266 + t100 - t54 + t94 + t95) - t217 
     #* t387 + t186 * t21 * t394) + t370 * (t201 * t143 + t197 + t8 * (t
     #21 * (t48 * t32 - t143) - t175) + t93 * (t170 * (t22 * (-t54 + t8)
     # - t4) + t144) + (t20 * (t102 + t55 + t57 + t101) + t130 + t217 * 
     #(-t100 - t94 - t95 - t55 + t57)) * t22 * t1)
      t197 = t146 * t14
      t209 = t14 * t20
      t230 = t17 - t1
      t305 = t151 * t17
      t326 = t101 * t20
      t329 = t305 * t121
      t127 = t380 * (t15 * (t14 * (-t102 * t151 * t121 + t104 * t186) + 
     #t4 * (t21 * (-t48 + t47) - t326) - t183 - t104 * t8 * t71) + t185 
     #+ t100 * (-t175 - t329 + t144)) + t380 * (-t305 * t178 * t121 - t9
     #5 * t175 + t186 * (t20 * (t22 * (-t51 * t17 + t344) - t56) + t89 *
     # (t1 * t42 * t22 - t41) - t190 * (t56 + t54)) + t266 * (t20 * (t27
     #1 - t56 + t100 + t101) - t257 * t71 + t272 * t230 * t42 + t102 * t
     #191 + t94 * t191 + t95 * t191) + t95 * t144) + t32 * t278 * (t127 
     #* t161 * t14 - t168 + t151 * t114 * (t14 * (t121 * (-t101 - t55 + 
     #t57) - t155) - t136) + t170 * (t186 * (t57 - t14 + t55) - t47 - t4
     #8 + t215 * (-t102 - t55 - t57) + t375 * t220) + t144 * t15 * (t102
     # + t14 + t55 - t57 + t101) + t172 * (t190 * (t101 + t100) - t209))
      t155 = t267 * t276
      t157 = t157 * t1
      t135 = t280 * (t20 * (t287 * t14 - t57 * t36) + t8 * (t32 * t218 *
     # t263 - t157 * t14 - t140) + t173 * t14 * (t20 * t230 + t373) + t1
     #30 * t4 * t230 + t100 * t107 * t230 + t102 * t107 * t230 + t94 * t
     #107 * t230 + t101 * t107 * t230 + t95 * t107 * t230) + t374 * (t20
     # * (t47 * t230 + t135) - t159 * (t112 + t382) + t198 * (t71 * (t14
     #8 - t110) - t340) + t215 * (t20 * (t17 * (t57 + t55) - t254) + t26
     #2 * (t201 * t15 * t372 + t151 * t230) + t410 * t372) + t364 * t230
     # + t288 * (t71 * t253 - t47 - t48))
      t148 = t376 * t42
      t183 = t17 * t121
      t115 = -t262 * t56 * t107 * t115 + t170 * t153 * t121 * t372 + t35
     #8 * t170 * t121 * t372 + t256 + t174 * t14 * (t175 + t329 - t144) 
     #+ t172 * t14 * (t175 - t329 + t144) + t279 * t146 - t209 * t136 * 
     #t115 + t143 * t15 * (t347 + t326) + t197 * t237 * (-t102 + t54) + 
     #t277 * t95 * t121 * (t14 * (t104 * t114 - t170) + t146)
      t115 = t263 * t115 + t212 * t263 * (-t340 * t93 + t71 * (t21 * (-t
     #48 * t8 + t36 * (t102 + t94 + t55 - t57 + t100 + t101)) + t93 * (t
     #8 * t270 + t410)) + t333 * (-t102 - t94 - t95 - t55 - t100 - t101)
     # + t157 * (t93 * t55 + t8 * (t100 + t94 + t95 - t54 + t55))) + t37
     #0 * (t137 * (t143 * t50 - t183 * t192 + t215 * t46) + t20 * (t183 
     #* (t42 * (t375 - t28) - t101) * t1 + t215 * (-t32 * t263 + t48) + 
     #t143 * (t102 + t148)) + t259 * (t20 * (t102 + t148 + t101) + t130)
     # + t186 * (t1 * (-t55 * t17 * t21 * t121 + t258) - t140) + t179 * 
     #t15 * t121 * (-t14 * t41 + t389))
      t115 = -t291 * t205 * t14 * t276 * t236 + t127 * t289 + t16 * t135
     # + t233 * t113 * t263 * t114 * t371 - t285 * t115 - t290 * (t155 *
     # (t266 * (t21 * (t56 + t94 + t95 - t55) + t201) - t209 + t21 * t15
     # * (-t4 * (t102 + t94 + t95 + t100) - t56 * t173)) + t155 * (t20 *
     # (t102 + t94 + t95 - t54 + t100) + t42 * (t15 * (-t262 * t132 + t3
     #52) + t337) + t130 + t186 * (t21 * (t57 - t14 - t54) + t201 + t190
     # * t210) + t55 * t221 + t101 * t221 + t388 * (t102 + t14 + t54 + t
     #100 + t101))) + t35 * t255 * t17 * t371 - t113 * (t8 * t71 * t372 
     #+ t204 + t262 * (t56 * (t114 * t14 - t277) + t119) + t48 * t372 + 
     #t47 * t372 + t93 * t71 * t372) + t146 * t237 * t245 * t276
      t127 = t265 * t25
      t135 = t95 * t265
      t140 = t100 * t265
      t143 = t93 * t265
      t148 = t102 * t265
      t153 = t101 * t265
      t155 = t94 * t265
      t157 = t146 * (-t102 - t10 - t55)
      t173 = -t368 - t303
      t175 = t146 * (-t57 - t8 - t95 - t93 - t101 - t100 - t94)
      t51 = t52 + t13 + t51
      t165 = t21 * (t6 * (-t102 - t93 - t94 - t95 + t8 + t57 - t100 - t1
     #01) - t344) - t165 - t264 * t366
      t179 = -0.1q1 / t79
      t183 = -0.1q1 / t31
      t79 = 0.1q1 / t79
      t27 = 0.1q1 / t27
      t31 = 0.1q1 / t31
      t191 = 0.1q1 / t39
      t26 = 0.1q1 / t26
      t92 = 0.1q1 / t92
      t198 = 0.1q1 / t45
      t38 = 0.1q1 / t38
      t201 = 0.1q1 / t28
      t204 = 0.1q1 / t97
      t205 = 0.1q1 / t81
      t37 = 0.1q1 / t37
      t96 = 0.1q1 / t96
      t33 = 0.1q1 / t33
      t40 = 0.1q1 / t40
      t209 = t2 * t26
      t210 = t92 * t39
      t215 = t210 + t209
      t220 = t209 * t191
      t230 = t220 + t92
      t233 = t204 ** 2
      t253 = t41 ** 2
      t254 = t90 ** 2
      t256 = t43 ** 2
      t257 = t78 * t37
      t258 = t257 * t76
      t259 = t258 * t79
      t205 = t259 * t205
      t264 = t23 * t38
      t265 = t264 * t19
      t266 = t265 * t31
      t198 = t266 * t198
      t267 = t2 * t24
      t270 = t158 * t40
      t271 = t40 * t215
      t272 = t32 * t40
      t277 = t27 * t3
      t74 = 0.1q1 / t74
      t83 = 0.1q1 / t83
      t72 = 0.1q1 / t72
      t73 = 0.1q1 / t73
      t62 = 0.1q1 / t62
      t70 = 0.1q1 / t70
      t87 = 0.1q1 / t87
      t279 = t65 * t68 * t73
      t280 = t86 * t72
      t285 = t280 * t58
      t287 = t285 - t279
      t288 = t37 ** 2
      t289 = t38 ** 2
      t290 = t34 ** 2
      t291 = t23 ** 2
      t305 = t80 ** 2
      t326 = t285 * t85
      t329 = t279 * t64
      t333 = t49 * t74
      t19 = t19 * t291
      t76 = t76 * t78 ** 2
      t136 = t119 * t3 * t136
      t339 = t67 * t82
      t344 = t339 * t333
      t287 = t277 * (t344 * t92 * t96 * t234 * t287 + (t344 * t2 * t96 *
     # t234 * t287 + t272 * (t136 * t29 * t33 * (-t19 * t31 * t289 * t29
     #0 + t76 * t79 * t305 * t288) + (t36 * (t258 * t359 * t80 * t70 * t
     #87 - t265 * t360 * t83 * t62 * t34) * t33 + t333 * (-t326 * t88 * 
     #t83 * t87 + t329 * t66 * t62 * t70)) * t234 * t67 * t2)) * t26 * t
     #191)
      t344 = 0.1q1 / t60
      t347 = 0.1q1 / t59
      t339 = t339 + t267
      t352 = t85 * t40
      t358 = t352 * t39 + t86
      t368 = t64 * t40
      t372 = t368 * t39
      t373 = t372 + t65
      t374 = t65 * t191
      t375 = t86 * t191
      t376 = t72 * t58
      t258 = t258 * t87 * t70
      t380 = t237 * t36
      t382 = t380 * t40
      t358 = t27 * (t333 * ((-t374 * t272 * t150 * t60 * t64 * t275 * t2
     #6 * t62 * t70 * t339 + t2 * (-t210 * t373 + t209 * (-t372 - t65)) 
     #* t96 * t344 * t63) * t73 * t68 + t376 * (t375 * t272 * t150 * t59
     # * t85 * t26 * t275 * t87 * t83 * t339 + t2 * (t209 * t358 + t210 
     #* t358) * t96 * t347 * t84)) + t382 * t33 * t26 * t234 * t191 * t1
     #50 * (t265 * t339 * t62 * t83 - t258 * t339))
      t372 = t21 * t91 - t90 * (-t132 * t35 + t89)
      t388 = -t210 - t209
      t389 = -t220 - t92
      t390 = t201 ** 2
      t391 = t183 ** 2
      t392 = t179 ** 2
      t393 = t388 * t40
      t394 = t237 * t3
      t397 = t394 * t392 * t391 * t191 * t26
      t398 = t43 * t3
      t399 = t398 * t29
      t400 = t220 * t20
      t404 = t5 * t67
      t405 = t67 * t4
      t406 = t333 * t40 * t96 * t3
      t359 = t406 * (-t404 * t210 * t376 + t209 * (t376 * t24 * t85 + (-
     #t24 * t64 + t404) * t73 * t68)) * t234 + t405 * (t399 * t205 * t28
     #6 * t359 * t26 * t80 * t33 * t191 * t40 + (t393 * t1 - t400) * t96
     # * t347 * t344 * t53)
      t407 = t2 * t5
      t408 = t64 * t82
      t409 = t376 * t85
      t410 = t44 * t90
      t82 = -t410 * t91 * t3 * t390 * t92 * t96 * t204 + t406 * (t209 * 
     #t67 * (-t408 * t68 * t73 + t376 * (t82 * t85 - t407)) + t210 * (t4
     #09 * t339 + (t67 * (-t408 + t407) - t267 * t64) * t73 * t68)) * t2
     #34 + t405 * (-t399 * t198 * t286 * t360 * t191 * t26 * t34 * t33 *
     # t40 - t20 * t53 * t92 * t344 * t347 * t96) * t2
      t82 = t27 * t82 + t27 * (t2 * t359 + t41 * (t399 * (t397 * t90 * t
     #281 * t29 * t40 + (t389 * t91 * t90 - t271 * t313) * t96 * t201) *
     # t233 + t398 * (-t397 * t235 * t40 * (t21 * t281 + t90 * (-t35 * (
     #t21 * (-t319 * t235 * t249 + t1 * (t10 + t8) * t29 * t241) - t255 
     #* t25) - 6 * t112 * t6 * (t107 - t268) - 3 * t28 * (t104 * t301 * 
     #t241 * t235 + t142 * t228 - t146 * t365) + 4 * t282 * t6 * (t21 * 
     #(-t13 * t241 * t235 + t36) + t354) + t28 * (-t104 * t119 + t93 * (
     #t30 * t367 + t146 + t300)) + t282 * t379 * (-8 * t1 + 5 * t142))) 
     #+ (-t220 * t90 * t91 + t393 * t313) * t96 * t390 + t29 * (t220 * t
     #372 + t372 * t92) * t96 * t201) * t204))
      t91 = -t42 * t86 + t84 * (t280 - t21)
      t107 = -t21 * t63 - t42 * t65
      t112 = t87 ** 2
      t142 = t33 ** 2
      t241 = t33 * t142
      t249 = t62 ** 2
      t255 = t62 * t249
      t268 = t73 ** 2
      t280 = t70 ** 2
      t282 = t74 ** 2
      t286 = t83 ** 2
      t300 = t23 * t18
      t313 = t300 * t115
      t339 = t37 * t112 * t280 * t80
      t114 = t339 * (t123 * (-t42 * t78 + t77 * (t257 - t21)) - t78 * t7
     #7 * (-t16 * t216 + t35 * t149 + t247 * (t187 * (t144 * t177 - t184
     # + t185 + (t20 * (-t57 + t54 - t55) - t130) * t15 * t4 + t100 * t1
     #76 + t102 * t176 + t101 * t176 + t95 * t176 + t186 * (t181 + t180)
     # + t151 * t8 * t15 * t121 * t178) + t116 * (-t169 * t114 - t170 * 
     #t164 - t168 + t172 * (t20 * (t102 + t94 + t95 + t54 - t56 + t100 +
     # t101) + t151 * (t102 + t55 - t57 + t100 + t101) * t22) + t174 * (
     #(t163 * t25 - t166) * t22 * t1 + t134 + t167))) - t246 * t206 + 3 
     #* t22 * t116 * t1 * (t1 * t218 + t199 * t221 + t222 - t223) + t245
     # * t238 * t21 * t236 + t244 * t212 * t117 * (t243 * t32 + t20) - t
     #145 + 4 * t113 * (t239 * t75 + t240) + t252 * t1 * (t15 * (t11 * (
     #-(0.5q1 / 0.2q1) * t128 + t248) - t250) + t251)))
      t116 = t65 * t63
      t8 = t116 * (t16 * t330 - t244 * t292 * t7 * t301 * (t20 * (-t52 -
     # t13) + t202) + t245 * t292 * t298 * t384 - t246 * t348 - t247 * t
     #182 - t35 * t229 - t301 * (t8 * t106 - t132 * t311 + t93 * t299 + 
     #t310) * t293 * t294 - t126 - (t20 * (t4 * (-t139 + t295) + t232 * 
     #t55) + t306 - t317 * t104) * t234 * t60 - t119 * t137 * t328 - t29
     #8 * (t10 * (t102 + t100 + t101) - t124) * t293 * t294 - (t119 * (t
     #21 * (-t100 * t4 + t28 * t311) + t105 * (t102 + t57) + t95 * t299)
     # + t297 + t231 * (t20 * (t102 + t57 - t54 + t101 + t100 + t94) + t
     #130)) * t234 * t60 - (t20 * (t120 + t122 * (-t54 + t55 + t100) + t
     #383 * t355) + t316 + t396 + t401 * (t55 + t95)) * t234 * t60 - t40
     #2 - t298 * (-t302 - t161 + t10 * (-t57 + t95 + t55 + t94) + t125) 
     #* t293 * t294 + 3 * t332 * t395 * t341 + 4 * t353 + t361)
      t13 = t86 * t84
      t6 = t13 * (-t16 * t211 + t244 * t378 * t7 * t301 * (t242 + t377) 
     #+ t245 * t378 * t298 * t384 + t246 * t362 + t247 * t208 + t35 * t3
     #35 - t1 * (t137 * (t46 * (t10 + t93) + t36 * t50) + t110 * t224 + 
     #t146 * (t4 * (t102 + t57 + t95 + t101 + t100 + t94) + t309)) - t29
     #8 * (-t10 * t103 - t124 + t125) * t293 * t296 - (t131 * t231 + t29
     #7) * t234 * t59 - t47 * t1 * t261 - (t307 * t299 + t48 * (t304 + t
     #303) + t146 * (t261 * t50 - t273 * t7 + t302) + t308 * (-t138 + t1
     #92) + t306) * t234 * t59 - t160 * t109 - t6 * (-t132 * t311 * t6 +
     # t106 * t319 + t299 * t320 + t312 * t321 - t322) * t293 * t296 - (
     #t20 * (t20 * (t315 + t314) + t161 * t1) + t316 + t57 * t119 * t106
     # + t317 * t299 + t318 * t299) * t234 * t59 - 4 * t325 * (t323 * t3
     #24 - t36) + 3 * t343 * t1 * t341 + t334)
      t7 = t36 * (t114 + (t313 * t34 * t289 + (t115 * (-t18 * t21 - t23 
     #* t42) - t300 * (t245 * t363 * t276 * t236 - t16 * t171 - t244 * t
     #212 * t276 * (t32 * t203 - t20) + t246 * t193 - t247 * (t263 * (-t
     #144 * t385 * t15 + t168 + t274 + t172 * (t190 * t214 - t20 * t387)
     # + t174 * (t20 * (-t57 - t55 + t56) + t130 + (t147 * t42 + t207) *
     # t22 * t1)) + t370 * (t101 * t260 + t184 + t340 + t228 * t32 * (-t
     #100 - t94 - t95 - t54 + t56) + t15 * (t20 * (t100 + t94 + t95 - t5
     #4 + t55 + t57) + t130) * t4 + t364 + t151 * (t1 * t381 + t54 * t19
     #9) * t121 + t102 * t260)) + t35 * t194 - t111 - 3 * t22 * t263 * t
     #1 * (t1 * (-t217 - t20) - t222 + t223 + t199 * t269) - 4 * t113 * 
     #(t17 * t323 - t240) + t197 * (t15 * (t14 * (-(0.5q1 / 0.2q1) * t26
     #2 + t248) - t250) + t251))) * t34 * t38) * t249 * t286)
      t10 = t282 * t49
      t15 = t237 * t40 * t26 * t275 * t191 * t27 * t150
      t22 = t257 * t77 * t123 * t80 * t280 * t112
      t22 = t15 * (t10 * (t285 * t112 * t286 * t188 * t84 * (t74 + t83 +
     # t87) + t279 * (t249 * (t70 * t280 - t280 * t74) + t255 * t280) * 
     #t108 * t63) + t36 * (t142 * (-t22 * (t70 + t87) + t313 * (-t83 * t
     #286 * t249 * t38 * t34 - t286 * t255 * t38 * t34)) + t241 * (t264 
     #* t18 * t115 * t286 * t249 * t34 - t22)))
      t47 = t272 * t277 * t209 * t234 * (t36 * (t265 * t18 * t83 * t62 *
     # t34 - t258 * t77 * t80) * t33 + t333 * (t326 * t84 * t83 * t87 - 
     #t329 * t63 * t62 * t70))
      t48 = t72 ** 2
      t52 = t344 ** 2
      t55 = t49 ** 2
      t57 = t53 ** 2
      t94 = t65 ** 2
      t95 = t31 ** 2
      t100 = t86 ** 2
      t101 = t79 ** 2
      t102 = t68 ** 2
      t103 = t347 ** 2
      t58 = t58 ** 2
      t105 = t107 * t219
      t106 = t116 * (-t351 * t141 * t21 + t35 * t20 * (t351 * t366 + t56
     # + t9) - t154 - t156 - t196)
      t109 = t376 * t347
      t110 = t109 * (t226 * t91 - t13 * (-t324 * t141 * t21 + t35 * t20 
     #* (t324 * t366 + t56 + t9) - t154 - t156 - t196))
      t111 = t264 * t284 * t34
      t113 = t257 * t129 * t80
      t114 = t272 * t3
      t115 = t382 * t3 * t235 * t26
      t117 = t115 * t142 * t191
      t119 = t10 * t96
      t72 = t61 * (t72 * t84 - t42) - t84 * (-t247 * t342 * t165 - t35 *
     # t146 * (t324 * t51 + t54 + t56 + t9) - t324 * (t41 * (t127 + t403
     #) - t135 - t140 - t143 - t148 - t153 - t155) - t157 - t324 * t173 
     #- t175)
      t42 = t42 * t189
      t51 = t63 * (-t247 * t332 * t165 - t35 * t146 * (t51 * t351 + t54 
     #+ t56 + t9) + t351 * (t41 * (-t127 - t403) + t135 + t140 + t143 + 
     #t148 + t153 + t155) - t157 - t351 * t173 - t175)
      t54 = t271 * t189
      t56 = t65 * t219
      t18 = t277 * (-t159 * t356 * t57 * t92 * t96 * t103 * t52 + t136 *
     # t272 * t26 * t33 * t191 * t234 * (t19 * t18 * t290 * t289 * t83 *
     # t62 - t76 * t77 * t305 * t288 * t87 * t70) + t119 * t234 * (t68 *
     # ((t92 * ((-t51 - t42) * t40 * t39 + t105 - t106) - t209 * t40 * (
     #t51 + t42)) * t344 * t73 + t63 * (t56 * t92 + t54) * t344 * t268) 
     #+ t109 * t40 * (t209 * t72 + t210 * t72))) + t277 * ((t393 * t146 
     #* t386 - t400 * t356) * t96 * t103 * t52 * t57 * t41 + t114 * t191
     # * t26 * (t85 * t100 * t84 * t58 * t83 * t48 * t87 - t64 * t94 * t
     #63 * t102 * t268 * t62 * t70) * t234 * t74 * t55 + t117 * (t113 * 
     #t79 * t101 + t111 * t31 * t95) + t119 * (t110 * t92 + t220 * (t68 
     #* (t116 * t219 * t344 * t268 + (-t106 + t105) * t344 * t73) + t110
     #)) * t234)
      t19 = t96 * t41
      t31 = t66 * t68 * t73 * t344
      t42 = t272 * t150 * t191
      t51 = t96 * t99
      t57 = t40 * t27 * t2 * (t404 * t380 * t264 * t150 * t26 * t33 * t1
     #91 * t234 * t83 * t62 + t333 * (t210 * t67 * t96 * (t409 * t88 * t
     #347 - t31 * t64) + t404 * t32 * t191 * t26 * (-t279 * t60 * t62 * 
     #t70 + t285 * t59 * t83 * t87) * t275 * t150 + t51 * (t279 * t215 -
     # t285 * t215) * t234 * t3))
      t5 = t57 + t27 * (t3 * (t394 * t40 * t191 * (t410 * t281 * t392 * 
     #t391 * t26 * t204 * t235 * (t179 + t183) + t36 * (-t291 * t83 * t6
     #2 * t38 + t257 * t70 * t87 * (-t209 * t5 + t78)) * t33 * t234 * t6
     #7) + t19 * (t389 * t347 * t344 * t53 * t137 + t271 * (t20 * (t152 
     #- t227) - t89 * t158) * t29 * t201 * t204 * t43 - t271 * t146 * t3
     #47 * t344 * t53)) + t333 * t67 * (t26 * (-t31 * (t374 + t368) + t1
     #09 * (t375 + t352) * t88) * t96 * t2 ** 2 + t92 * (-t279 * t66 * t
     #344 + t285 * t88 * t347) * t96 * t2 + t42 * t275 * (t60 * t94 * t6
     #8 * t73 * t62 * t70 - t376 * t59 * t100 * t83 * t87)))
      t31 = t86 * t226
      t49 = t277 * t96 * t234 * t74 * t282 * t49 * (t109 * t84 * (t31 * 
     #t230 + t271 * t61) + (t56 * t389 - t54) * t344 * t73 * t68 * t63)
      t45 = -0.1q1 / t45
      t54 = 0.1q1 / t41
      t57 = -0.1q1 / t81
      t59 = 0.1q1 / t97 ** 2
      t51 = t51 * t388
      t60 = t32 * t67
      t62 = t96 * t230 * t24
      t54 = t96 * t54 * t67
      t41 = t256 * t41
      t64 = t270 * t32
      t24 = t27 * t2 * (t54 * t201 * t20 * (-t230 * t4 + t271 * t46) + t
     #41 * (-t64 * t150 * t24 * t28 * t179 * t183 * t191 * t26 * t235 + 
     #(t40 * (t60 * t158 * t331 * t179 * t183 * t191 * t26 + t51) + t62)
     # * t29 * t3 + t96 * t201 * t67 * (-t230 * t331 + t271 * (t30 * t50
     # + t46))) * t233 * t90 + (-t41 * t42 * t24 * t98 * t45 * t26 * t57
     # + t398 * (t40 * (t60 * t98 * t327 * t45 * t191 * t26 * t57 + t51)
     # + t62) + t54 * (-t230 * t327 + t271 * (t44 * t50 - t46))) * t59 *
     # t235 * (t21 * t44 + t20) * t28)
      t30 = t1 * t39 * t40
      t41 = t93 * t92
      t42 = t267 * t230
      t2 = t27 * (t2 * (t19 * (-t210 * (t30 + t20) + t209 * (-t30 - t20)
     #) * t347 * t344 * t53 + t114 * t36 * t29 * t26 * t33 * (-t259 * t8
     #0 + t266 * t34)) + t333 * t96 * t3 * (t285 * ((-t41 - t2) * t40 * 
     #t67 - t42) + t279 * ((t41 + t2) * t40 * t67 + t42)) * t234)
      t4 = t117 * (t34 * t38 * t95 * (t23 * (-t16 * t200 * t278 * t162 +
     # t246 * t370 * t350 - t247 * t370 * t357 + t35 * t1 * (t28 * (t20 
     #* (t9 + t14) + t104 * t17) + t132 * (t262 * t349 + t36)) + t20 * (
     #-t152 * t17 * t4 + t369 * t278 - t336) - t345 * t371 - t346 * (t71
     # * t21 + t20) + t93 * (t17 * (-t133 + t337) - t338) - t195 * t17 *
     # (3 * t262 + t12)) + t284 * (-t264 + t21)) + (t37 * (t21 * t129 + 
     #t78 * (-t16 * t200 * t118 * t162 + t246 * t187 * t350 - t247 * t18
     #7 * t357 + t35 * t1 * (t28 * (t11 * t218 + t354) + t132 * (t128 * 
     #t349 + t36)) - t20 * (t283 * t152 - t213 * t25 + t336) - t345 * t2
     #25 - t346 * (t69 * t21 + t20) - t93 * (t75 * (t133 - t337) + t338)
     # - t195 * t75 * (3 * t128 + t12))) - t78 * t129 * t288) * t80 * t1
     #01)
      t4 = t277 * (t4 + (t230 * t347 * t48 * t58 * t84 * t100 + t13 * t2
     #71 * t85 * t347 * t48 * t58 + t116 * t344 * t268 * t102 * (-t373 *
     # t92 + t209 * (-t374 - t368))) * t96 * t74 * t55)
      t9 = t61 * t40
      t11 = t189 * t40
      t12 = t376 * t84 * t103
      t13 = t68 * t63 * t73 * t52
      t9 = t277 * (t115 * t241 * t191 * (-t113 * t101 + t111 * t95) + t1
     #19 * (t92 * (t13 * (t11 * t39 + t56) + t12 * (t9 * t39 + t31)) + t
     #209 * (t12 * (t375 * t226 + t9) + t13 * (t374 * t219 + t11))))
      ret = -64 * t287 + 96 * t358 - 16 * t82 + 2048 * t15 * (t7 * t142 
     #+ t10 * (t68 * (t116 * t108 * t249 * t280 * t268 + (t107 * t108 - 
     #t8) * t280 * t249 * t73) + t376 * t112 * t286 * (t188 * t91 - t6))
     #) + 8192 * t22 - 320 * t47 - 256 * t18 - 32 * t5 - 1024 * t49 + 48
     # * t2 + 128 * t4 - 512 * t9 - 8 * t277 * (t272 * t43 * t29 * t26 *
     # t191 * t3 * (-t158 * t254 * t253 * t179 * t183 * t233 * t43 + t26
     #7 * t36 * t32 * t33 * (t198 - t205)) + (t20 * (t271 * t1 + t20 * t
     #230) + t253 * (t270 * t215 * t256 * t233 * t90 + t230 * t256 * t23
     #3 * t254)) * t96 * t201) + 4 * t24 + 12 * t64 * t410 * t277 * t209
     # * t29 * t179 * t183 * t204

      hjetmass_qqbgg_bubble_pmmp_s123 = -ret/16q0*(0,1q0)
      return

      end function

      complex*32 function hjetmass_qqbgg_bubble_pmmp_s124
     &     (i1,i2,i3,i4,za,zb,mt,p,flip)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt
          real*16 p(mxpart,4)
          complex*32 alpha
          logical flip

      alpha = (za(i1,i2)*zb(i2,i1) + za(i1,i4)*zb(i4,i1) +
     & za(i2,i4)*zb(i4,i2))/(za(i1,i2)*zb(i2,i1) + za(i1,i4)*zb(i4,i1))

      p(5,:) = real(alpha)*p(i1,:)
      p(6,:) = (1q0-real(alpha))*p(i1,:) + p(i2,:) + p(i4,:)
      if (flip .eqv. .true.) then
          call spinoru_qp(6,p,zb,za)
      else
          call spinoru_qp(6,p,za,zb)
      end if

      t1 = za(i2, i3)
      t2 = zb(i1, 6)
      t3 = za(6, i2)
      t4 = zb(i2, 5)
      t5 = za(6, i3)
      t6 = zb(i3, 5)
      t7 = za(6, i4)
      t8 = zb(i4, 5)
      t9 = t3 * t4
      t10 = t5 * t6
      t11 = t7 * t8
      t12 = t11 + t9 + t10
      t13 = za(5, i2)
      t14 = za(6, i1)
      t15 = zb(i2, 6)
      t16 = za(5, i3)
      t17 = zb(i3, 6)
      t18 = za(5, i4)
      t19 = zb(i4, 6)
      t20 = t13 * t15
      t21 = t16 * t17
      t22 = t18 * t19
      t23 = t22 + t20 + t21
      t24 = t14 * t2
      t25 = t13 * t4
      t26 = t3 * t15
      t27 = t16 * t6
      t28 = t5 * t17
      t29 = t18 * t8
      t30 = t7 * t19
      t31 = t24 - t25 + t26 - t29 + t30 - t27 + t28
      t31 = t31 ** 2
      t32 = 4 * t12 * t23 + t31
      t32 = sqrt(t32)
      t33 = t24 - t25 + t26 - t29 + t30 - t27 + t28 + t32
      t34 = 0.1q1 / t23
      t35 = (0.1q1 / 0.2q1)
      t36 = t35 * t13
      t37 = t36 * t33 * t34 - t3
      t38 = t35 * t15
      t39 = t38 * t33 * t34 + t4
      t40 = zb(i4, i1)
      t41 = za(i1, i2)
      t42 = za(i2, i4)
      t43 = zb(i2, i1)
      t44 = zb(i3, i1)
      t31 = 4 * t12 * t23 + t31
      t31 = sqrt(t31)
      t45 = t24 - t25 + t26 + t31 - t29 + t30 - t27 + t28
      t46 = 0.1q1 / t18
      t47 = t35 * t33 * t34
      t48 = t7 * t46
      t49 = -t48 + t47
      t50 = 0.1q1 / t19
      t51 = t8 * t50
      t52 = t51 + t47
      t45 = 0.1q1 / t45
      t53 = 2 * t12
      t54 = t53 * t45 + t47
      t32 = t24 - t25 + t26 - t29 + t30 - t27 + t28 - t32
      t55 = t34 * (-t32 + t33)
      t36 = t36 * t32 * t34 - t3
      t56 = zb(i4, i2)
      t31 = t24 - t25 + t26 - t31 - t29 + t30 - t27 + t28
      t57 = t35 * t32 * t34
      t58 = -t48 + t57
      t59 = t51 + t57
      t31 = 0.1q1 / t31
      t53 = t53 * t31 + t57
      t38 = t38 * t32 * t34 + t4
      t60 = t57 * t16 - t5
      t61 = t47 * t16 - t5
      t62 = t13 * t7
      t63 = t62 * t46
      t64 = t63 - t3
      t65 = t7 * t16
      t66 = t65 * t46 - t5
      t67 = t48 + t51
      t68 = za(i1, i4)
      t69 = t41 * t43
      t70 = t68 * t40
      t71 = t42 * t56
      t72 = t71 + t69 + t70
      t73 = t26 + t24
      t74 = t3 * t17
      t75 = t6 * t13
      t76 = -t75 + t74
      t77 = t10 + t9
      t78 = t13 ** 2
      t79 = t46 ** 2
      t80 = t7 ** 2
      t81 = t3 ** 2
      t82 = t16 ** 2
      t83 = t24 * t16
      t84 = t6 * t82
      t85 = t81 * t15
      t86 = t85 * t16
      t87 = t13 * t5
      t88 = t5 * t78
      t89 = t88 * t15
      t81 = t81 * t4
      t90 = t81 * t16
      t91 = t13 * t16
      t92 = t91 * (t21 + t20)
      t93 = t5 * t3
      t94 = t27 + t25
      t95 = -t10 - t9
      t96 = t93 * t17
      t97 = t46 * t80
      t76 = t7 * (-t92 * t79 ** 2 * t80 ** 2 + t93 * t77 + t48 * (t5 * (
     #t3 * t73 + t5 * t76) - t90) + (t16 * (t13 * (-t25 + t24) + t16 * t
     #76) + t89) * t46 * t79 * t7 * t80 + (t3 * (t84 - t83) - t86 + t87 
     #* (-t28 - t24 + t25)) * t79 * t80) - 2 * t97 * (t48 * (t16 * (t13 
     #* t95 + t96) + t20 * t93) - t91 * (t28 + t26) * t79 * t80 + t93 * 
     #t94)
      t98 = za(i3, i4)
      t99 = za(i1, i3)
      t100 = 0.1q1 / t23
      t101 = t16 * t3
      t102 = t87 + t101
      t103 = t62 * t16
      t104 = t18 * t102
      t105 = t104 + t103
      t106 = t33 ** 2
      t107 = t33 * t106
      t108 = t100 ** 2
      t109 = t108 ** 2
      t110 = t100 * t108
      t111 = t18 * t3
      t112 = t99 * t44
      t113 = t98 * zb(i4, i3)
      t114 = t1 * zb(i3, i2)
      t115 = t110 * t107
      t116 = t84 * t7
      t117 = t5 * t8
      t118 = t18 * t5
      t119 = t111 * t16
      t120 = t119 * (t22 + t21)
      t121 = t74 * t82
      t122 = t11 + t9 + t10
      t123 = -t30 + t27
      t124 = t65 * t6
      t125 = t62 * t15
      t126 = t22 + t20
      t127 = t118 * t126
      t128 = t82 * t7
      t129 = t128 * t17
      t130 = (t18 * (t26 + t113 + t114 + t112 + t69 + t70 + t71 - t24) -
     # t125) * t16
      t131 = t130 + t127 - t129
      t132 = t8 * t13
      t133 = t3 * t19
      t134 = t75 + t74
      t135 = t111 * t8
      t136 = t91 * t80
      t137 = t136 * t8
      t138 = t93 * t18
      t139 = t34 * t7 * t5 * (t3 * (t26 + t113 + t114 + t112 + t69 + t70
     # + t71 + t24 - t25) + t5 * t134 + t7 * (t133 + t132))
      t96 = t96 * t16
      t140 = t96 * t18
      t141 = (t138 * (-t28 - t24 - t26) - t137 + t65 * (t5 * (-t75 - t74
     #) + t135)) * t108
      t142 = t102 * t7 + t138
      t143 = t26 - t24
      t144 = t9 * t16
      t145 = t70 * t142
      t146 = t69 * t142
      t147 = t113 * t142
      t148 = t108 * t106
      t149 = t13 * t80
      t150 = t24 * t62
      t151 = t62 * t5
      t90 = t148 * (t3 * (t114 * t65 - t116) - t90 * t18 + t16 * (t149 *
     # t19 + t111 * (t29 + t27) + t150) * t100 * t33 + t151 * (-t30 + t1
     #14 - t28)) + t148 * (t138 * t114 + t145 + t146 + t147 + t65 * (t18
     # * (-t133 - t132) + t87 * t17) * t100 * t33 + t112 * t142 + t71 * 
     #t142 + t62 * (t143 * t5 + t144))
      t152 = t29 + t25
      t153 = t148 * t138
      t154 = t22 + t20
      t155 = (0.3q1 / 0.4q1)
      t156 = -(0.3q1 / 0.8q1)
      t157 = (0.1q1 / 0.4q1)
      t158 = (0.1q1 / 0.8q1)
      t159 = (0.1q1 / 0.16q2)
      t160 = 2 * t11 * t138
      t86 = -t35 * t33 * (t140 * t106 * t110 + t141 * t33 + t139) - t158
     # * (t115 * (t13 * (-t117 * t18 ** 2 + t116) + t69 * t104 - t120 * 
     #t100 * t33 + t121 * t7 + t86 * t18 + t4 * (-t118 + t65) * t78) + t
     #115 * (t16 * (t13 * (t7 * (t22 * t33 * t100 - t26 + t69) - t111 * 
     #t4) + t111 * t24) + t112 * t105 + t113 * t105 + t70 * t105 + t71 *
     # t105 + t114 * t105)) + t159 * t109 * t106 ** 2 * t13 * t131 + t15
     #7 * t90 - t155 * t153 * t152 + t156 * t138 * t115 * t154 - t93 * (
     #t123 * t18 * t108 * t106 + t122 * t7 + (t18 * t95 - t124) * t100 *
     # t33) + t160 * t33 * t100
      t90 = t48 * t19 + t8
      t95 = t5 ** 2
      t105 = t97 * t13
      t106 = t111 * t6
      t161 = t71 + t113 + t114 + t112 + t69 + t70
      t162 = t161 * t3
      t88 = t88 * t4
      t163 = t57 * t19 + t8
      t164 = t87 + t101
      t165 = t113 * t164
      t166 = t93 * (-t71 - t70 - t113)
      t167 = (t8 * (t18 * (-t87 - t101) + t103) + t70 * t164 + t71 * t16
     #4 + t69 * t164 + t114 * t164 + t112 * t164 + t165) * t100
      t84 = t84 * t3
      t168 = t93 * (-t24 - t112 - t26 - t30 - t28 - t69 - t114)
      t169 = (-t84 - t88) * t100
      t170 = t22 + t20 + t21
      t171 = t170 * t100
      t126 = t87 * t126
      t172 = (t13 * (t30 - t113 - t114 - t112 - t69 - t70 - t71 + t24) -
     # t111 * t19) * t16
      t173 = t172 - t121 - t126
      t174 = t32 ** 2
      t175 = t32 * t174
      t176 = t108 * t174
      t177 = -t155 * t176 * t173 - t167 * t32 - t169 * t32 - t166 - t168
     # - 2 * t93 * (t171 * t32 + t25 + t27 + t29)
      t34 = (t30 + t113 + t114 + t112 + t69 + t70 + t71 + t24 + t26 + t2
     #8) * t34
      t178 = (t29 + t25 + t27) * t100
      t179 = t71 * t102
      t180 = t87 * (t29 + t25)
      t121 = -t121 * t110 - t126 * t110 + t172 * t110
      t126 = -t157 * t176 * (t112 * t102 + t113 * t102 + t114 * t102 + t
     #69 * t102 + t70 * t102 + t179 - t180 - t84 + t8 * (t62 - t111) * t
     #16) - t158 * t175 * t121 + t35 * t93 * t32 * ((-t22 - t20 - t21) *
     # t108 * t32 + t34) - t93 * (t178 * t32 - t10 - t11 - t9)
      t75 = (t75 + t74) * t5
      t172 = t93 * (t18 * (t10 + t9) + t124)
      t135 = (t65 * (t75 - t135) + t138 * (t28 + t24 + t26) + t137) * t1
     #00
      t137 = -t87 - t101
      t181 = t137 * t18 - t103
      t182 = t85 * t18
      t78 = t4 * t7 * t78
      t183 = t180 * t18
      t128 = t128 * t134
      t134 = t114 * t181
      t184 = t71 * t181
      t185 = t70 * t181
      t186 = t113 * t181
      t187 = t25 - t24
      t188 = t69 * t181
      t181 = t112 * t181
      t189 = t111 * (-t30 + t27 + t29)
      t190 = t62 * (t30 + t24 + t28 - t29) + t189
      t191 = t176 * t138
      t127 = t127 * t110 - t129 * t110 + t130 * t110
      t129 = -t164 * t7 - t138
      t130 = t7 * t5
      t164 = t118 + t65
      t192 = t70 * t129
      t103 = t22 * t103
      t193 = t138 * t100
      t123 = t123 * t100
      t194 = -(0.9q1 / 0.4q1)
      t195 = (0.3q1 / 0.2q1)
      t57 = t57 * t18 - t7
      t196 = t47 * t18 - t7
      t47 = t47 * t19 + t8
      t197 = 0.1q1 / t17
      t198 = t6 * t197
      t199 = t198 * t19 - t8
      t200 = 0.1q1 / t16
      t201 = t5 * t200
      t202 = t198 + t201
      t203 = t198 * t13
      t204 = t3 + t203
      t205 = t87 * t200 - t3
      t206 = -t7 * t64
      t207 = t48 + t198
      t208 = -t62 + t111
      t34 = t157 * t148 * (t8 * t208 * t16 + t112 * t137 + t113 * t137 +
     # t114 * t137 + t69 * t137 + t70 * t137 + t71 * t137 + t180 + t84) 
     #- t158 * t107 * t121 - t35 * t93 * t33 * (t170 * t108 * t33 - t34)
     # - t93 * (t178 * t33 - t10 - t11 - t9)
      t121 = -t155 * t148 * t173 - 2 * t93 * (t171 * t33 + t25 + t27 + t
     #29) - t167 * t33 - t169 * t33 - t166 - t168
      t166 = t118 + t65
      t91 = t9 * t91
      t167 = t13 * t95 * t17
      t168 = t87 * t80 * t19
      t169 = t110 * t175
      t160 = t160 * t32 * t100
      t74 = -t155 * t191 * t152 + t156 * t169 * t138 * t154 + t157 * (t1
     #76 * (t7 * (t13 * (t83 * t32 * t100 + t112 * t5) + t179) + t147 + 
     #t114 * t142) + t176 * (t7 * (t87 * t143 - t167 - t84 + t91) + t145
     # + t146 - t168 + t111 * (t71 * t5 - t144) + t16 * (t62 * (t30 + t2
     #8 - t29) + t189) * t100 * t32 + t112 * t3 * t166)) + t158 * (t169 
     #* (-t69 * t87 * t18 + t186 + t101 * (t18 * (t25 - t69 - t24) + t12
     #5) - t71 * t13 * t166) + t169 * (t16 * (-t111 * t71 - t182 - t62 *
     # (t25 + t69)) - t128 + t134 + t181 + t183 + t185 + t18 * t16 * (t7
     #4 * t16 + t19 * t208) * t100 * t32)) + t159 * t109 * t174 ** 2 * t
     #13 * t131 - t35 * t32 * (t140 * t174 * t110 + t141 * t32 + t139) +
     # t93 * (t7 * (-t11 - t9 - t10) + (t18 * t77 + t124) * t100 * t32 +
     # t176 * t18 * (t30 - t27)) + t160
      t77 = t137 * t7 - t138
      t83 = t100 * t33
      t109 = -t48 + t201
      t124 = t132 * t50 + t3
      t131 = 0.1q1 / t67
      t137 = 0.1q1 / t6
      t72 = 0.1q1 / t72
      t139 = 0.1q1 / t33
      t140 = 0.1q1 / t32
      t141 = 0.1q1 / t41
      t55 = 0.1q1 / t55
      t44 = 0.1q1 / t44
      t53 = 0.1q1 / t53
      t142 = 0.1q1 / t207
      t54 = 0.1q1 / t54
      t143 = 0.1q1 / t202
      t145 = 0.1q1 / t43
      t68 = 0.1q1 / t68
      t146 = 0.1q1 / t42
      t147 = t40 * t141
      t43 = t43 * t68
      t152 = t147 * t145
      t156 = t152 + t68
      t157 = t40 * t156
      t158 = t157 * t98
      t159 = t1 * (t43 + t147) - t158
      t160 = t40 * t145
      t41 = t41 * t68
      t166 = t41 + t160
      t169 = t131 ** 2
      t170 = t64 * t206 * t145
      t171 = t93 * t72
      t173 = t171 * t143 * t200
      t174 = t37 * t45 * t54
      t178 = t36 * t53
      t180 = t178 * t31
      t189 = t72 * (t196 * t45 * t54 - t57 * t31 * t53) * t146
      t208 = t173 * t199
      t209 = t157 * t1
      t210 = -0.1q1 / t58
      t207 = 0.1q1 / t207
      t211 = 0.1q1 / t52
      t58 = 0.1q1 / t58
      t212 = -0.1q1 / t49
      t49 = 0.1q1 / t49
      t213 = 0.1q1 / t59
      t214 = t210 ** 2
      t215 = t36 ** 2
      t216 = t212 ** 2
      t217 = t37 ** 2
      t218 = t64 * t131
      t201 = t201 * t204
      t219 = t12 * t1
      t220 = t219 * t40 * t55
      t20 = t2 * (t201 * (2 * t106 * t197 + t7 * (t3 - t203)) * t207 ** 
     #2 * t143 * t197 + t51 * t108 * t216 * t214 * t131 * (-t64 * (3 * t
     #7 * (-t136 * t94 * t79 + t48 * (-t101 * t73 + t88) + t93 * (t28 + 
     #t24 + t26)) - 5 * t3 * t7 * (-t82 * t80 * t17 * t79 - t63 * t16 * 
     #t4 + t87 * t4) - 6 * t65 * t5 * (-t149 * t17 * t79 + t3 * t6) + 7 
     #* t97 * t20 * t3 * t66 + 2 * t5 * (t149 * t90 + t81 * t18 - t7 * (
     #t29 + t30) * t3) + 2 * t95 * (-t105 * t17 + t106) + 2 * t105 * t24
     # * t66 + 2 * t65 * (t7 * t90 * t3 - t105 * t90 - t81) - t7 * (t5 *
     # (t87 * t6 - t162) + t136 * (-t71 - t113 - t114 - t112 - t69 - t70
     #) * t79 + t48 * (t162 * t16 + t87 * t161)) + 4 * t97 * (-t92 * t79
     # * t80 + t27 * t102 + t48 * t89) - 8 * t96 * t97) + t76 * (t218 - 
     #t13))) * t146 * t79
      t27 = t44 * t2
      t9 = t27 * (t41 * t220 * t189 * t100 + (t20 + t220 * t50 * t100 * 
     #(-t217 * t39 * t211 * t45 * t54 * t49 + t215 * t38 * t58 * t53 * t
     #31 * t213) * t146 * t46 + t208 * t40 * t137) * t145 * t141) + t44 
     #* (t2 * (t208 * t68 * t137 + (t189 * t145 * t100 * t55 * t40 ** 2 
     #+ t152 * (t180 - t174) * t146 * t100 * t55) * t12 * t1 + t146 * (t
     #8 * (-t170 * t141 * t50 * t79 * t142 * t169 + (-t170 * t142 ** 2 +
     # t145 * (t13 * t206 - t64 * (-2 * t111 + t62)) * t142) * t79 * t50
     # * t141 * t131) - t173 * t166) * t197 * t2) + t72 * (t8 * (t159 * 
     #t3 + t209 * t7) + t157 * t9 * t1) * t140 * t139 * t23)
      t20 = -t41 - t160
      t63 = t160 * t98
      t65 = -t63 + t1
      t73 = t36 * t156
      t80 = t37 * t156
      t81 = t217 * t141
      t82 = t53 * t31
      t88 = t157 * t36
      t89 = t54 * t45
      t90 = t50 * t141
      t92 = t90 * t79
      t38 = t44 * (t8 * t2 * (t171 * t23 * t140 * t139 * t156 + t218 * (
     #-t92 * t216 * t108 * t145 * t146 * t210 * t214 - t92 * t212 * t216
     # * t108 * t145 * t146 * t214) * t76 * t2) + t72 * (t82 * (t140 * (
     #t1 * (t163 * (t36 * (-t43 - t147) - t157 * t57) - t88 * t38) + t15
     #8 * t36 * t163) + t88 * t100 * t99 * t2) + t89 * (t139 * (t47 * (t
     #159 * t37 + t209 * t196) + t209 * t37 * t39) - t157 * t37 * t100 *
     # t99 * t2)) * t55 * t12)
      t39 = t64 ** 2
      t43 = t143 ** 2
      t57 = t200 ** 2
      t64 = t66 * t210 * t212 * t100
      t66 = t142 * t197
      t76 = t95 * t204 * t57 * t43
      t76 = t27 * (t72 * t137 * (t76 * t156 * t199 - t3 * t8 * t156) + (
     #-t76 * t197 * t72 * t166 + t145 * (t95 * t204 ** 2 * t57 * t43 * t
     #207 * t197 + (-t66 + t64) * t50 ** 2 * t169 * t39 * t8 ** 2 + t219
     # * t50 * t108 * (-t89 * t33 * t217 * t211 * t49 + t82 * t32 * t215
     # * t58 * t213) * t55 * t56) * t46 * t141) * t146 * t2)
      t88 = 0.1q1 / t109
      t92 = 0.1q1 / t109
      t59 = -0.1q1 / t59
      t52 = -0.1q1 / t52
      t94 = 0.1q1 / t202 ** 2
      t67 = 0.1q1 / t67 ** 2
      t96 = t98 * t204
      t97 = t1 * (t198 * t18 + t7)
      t102 = t40 * t72
      t41 = t41 * t72
      t14 = t1 * t14 * t72 * t156
      t105 = t152 * t51
      t106 = t98 * t205
      t109 = t1 * (t118 * t200 - t7)
      t4 = t27 * (t197 * (-t105 * t39 * t98 * t200 * t46 * t88 * t146 * 
     #t142 * t131 + t43 * t40 * t5 * (t146 * (t145 * (t204 * t141 * (t96
     # * t207 * t46 + t1) + t102 * (-t96 - t97)) - t41 * (t96 + t97)) - 
     #t14 + t99 * t204 * t72 * t156) * t57) + t94 * t40 * t6 * (t146 * (
     #t145 * (t205 * t141 * (t106 * t92 * t46 + t1) + t102 * (-t106 - t1
     #09)) - t41 * (t106 + t109)) + t14 + t205 * t72 * t156 * t99) * t20
     #0 * t197 ** 2 + t90 * t146 * t145 * t100 * t46 * t1 * (t48 * (t40 
     #* (t51 * t15 - t4) - t51 * t2 * t56) * t52 * t67 * t59 * t124 ** 2
     # + t51 * t39 * t169 * t210 * t212 * (t40 * (t48 * t15 + t4) - t48 
     #* t2 * t56)))
      t6 = t27 * t46
      t14 = -t121 * t47 + t34 * (t47 * t54 - t19)
      t15 = t58 ** 2
      t43 = t140 ** 2
      t48 = t49 ** 2
      t51 = t139 ** 2
      t52 = t53 ** 2
      t57 = t55 ** 2
      t19 = t19 * t126
      t59 = t163 * t177
      t67 = t163 * t126
      t88 = t67 * t156
      t90 = t89 * t139
      t14 = t27 * (-t8 * t93 * t122 * t23 ** 2 * t43 * t51 * t72 * t156 
     #+ (t72 * (t140 * (t88 * t31 * t52 + t82 * (-t68 * (t59 + t19) + t1
     #52 * (-t59 - t19))) + t90 * (t152 * t14 + t14 * t68)) * t100 + t2 
     #* t141 * t79 * t145 * t146 * (t180 * t74 * t58 * t15 + t174 * t86 
     #* t49 * t48) * t108) * t57 * t12)
      t19 = t72 * t166
      t59 = -t37 * t86 * t141 * t79 * t145 * t48 + t19 * t34
      t43 = t27 * t57 * t12 * (t72 * (-t89 * t47 * t34 * t51 * t156 + t6
     #7 * t82 * (-t152 - t68) * t43) + (t89 * t59 + t82 * (t36 * t74 * t
     #15 * t141 * t79 * t145 + t20 * t72 * t126)) * t146 * t108 * t55 * 
     #t2)
      t51 = t54 ** 2
      t67 = t45 ** 2
      t31 = t31 ** 2
      t10 = t145 * (t36 * (2 * t135 * t32 + t155 * (t176 * (t101 * (t18 
     #* t187 + t125) + t188 + t181) + t176 * (t16 * (-t78 - t182) + t183
     # - t128 + t134 + t184 + t185 + t186)) + t35 * t175 * t13 * t127 + 
     #t195 * t176 * t16 * t190 + t194 * t191 * t154 + 2 * t172 - t130 * 
     #(t13 * (t10 + t11) + t3 * (t30 + t28) + t85) + t119 * (t22 + t21) 
     #* t110 * t175 - (t62 * t95 * t17 + t113 * t129 + t116 * t3 + t69 *
     # t129 + t151 * (-t71 - t26 + t30 - t114) + t144 * (-t62 + t111)) *
     # t100 * t32 - t103 * t175 * t110 - (-t114 * t164 * t3 - t71 * t3 *
     # t164 + t112 * t129 + t150 * t5 + t192) * t100 * t32 - t93 * t7 * 
     #(t24 + t71 + t112 - t25 + t69 + t70 + t113 + t114) - 3 * t193 * t3
     #2 * (t21 * t32 * t100 + t25 + t29) - 4 * t138 * (t123 * t32 - t11)
     #) + (-t178 + t13) * t74) * t79 * t141 * t15
      t10 = t82 * (t72 * (t166 * t53 * t126 - t166 * t177) + t10)
      t3 = (t45 * (t51 * t59 + t54 * (-t19 * t121 + t145 * (t13 * t86 + 
     #t37 * (t155 * (t148 * (t16 * (t111 * t187 - t70 * t62) + t184 + t1
     #86) + t148 * (t16 * (t125 * t3 - t182 - t78) - t128 + t134 + t181 
     #+ t183 + t188 - t104 * t70)) + t194 * t153 * t154 + t195 * t148 * 
     #t16 * t190 + t35 * t107 * t13 * t127 - t83 * (t7 * (t87 * (-t26 + 
     #t24) + t84 - t91 + t167) + t168 - t179 * t7 + t69 * t77 + t112 * t
     #77 + t114 * t77 - t148 * t120 + t111 * (-t113 * t5 + t144)) - t103
     # * t115 - t130 * (t3 * (t26 + t113 + t114 + t112 + t69 + t70 + t71
     # + t24 - t25) + t7 * (t133 + t132) + t75) - t83 * (-t138 * t71 - t
     #165 * t7 + t192) + 2 * t135 * t33 + 2 * t172 - 3 * t193 * t33 * (t
     #83 * t21 + t25 + t29) - 4 * t138 * (t123 * t33 - t11))) * t48 * t7
     #9 * t141)) + t10) * t146 * t108 * t55 * t2
      t5 = t72 * t12
      t7 = t27 * t55 * t12
      t10 = t90 * t156
      t11 = t5 * t27 * t100 * t55 * t57 * (t10 * t34 * t47 - t88 * t82 *
     # t140)
      t10 = t102 * t44 * t55 * t42 * t12 * (-t82 * t60 * t163 * t140 * t
     #156 + t10 * t47 * t61)
      t13 = t69 * t68 + t40
      t13 = t7 * t146 * t100 * (t1 * t72 * (-t174 * t13 + t180 * t13) + 
     #(t81 * t61 * t67 * t46 * t51 * t145 * t49 + t31 * t52 * t60 * t36 
     #* (t145 * (-t36 * t58 * t141 * t46 + t102) + t41) - t19 * t51 * t6
     #7 * t61 * t37) * t12 * t2)
      t15 = t6 * t152 * t12 * t100 * t55 * (-t174 * t61 * t49 + t180 * t
     #60 * t58)
      ret = 16 * t9 - 48 * t44 * (t157 * t117 * t42 * t23 * t140 * t139 
     #* t72 + ((t72 * (t40 * (t166 * t146 * t98 * t37 + t166 * t61) - t8
     #0 * t56 * t1) + t81 * t146 * t65 * t49 * t46) * t54 * t45 + t82 * 
     #(t72 * (t40 * (t20 * t146 * t98 * t36 - t166 * t60) + t73 * t56 * 
     #t1) + t215 * t141 * t146 * (t63 - t1) * t46 * t58)) * t100 * t55 *
     # t12 * t2) - 32 * t38 - 12 * t6 * t141 * (t8 * (t210 * t212 * t100
     # * t146 * t65 * t50 * t131 * t39 + t66 * t218 * t160 * t50) + t160
     # * t201 * t143 * t207 * t197) + 256 * t14 - 512 * t43 - 128 * t7 *
     # (t5 * (t80 * t51 * t67 * t139 * t47 * t61 - t73 * t60 * t163 * t3
     #1 * t140 * t52) + t3) + 1024 * t11 + 96 * t10 - 64 * t13 - 80 * t1
     #5 - 8 * t76 - 4 * t4 + 20 * t6 * t105 * t64 * t218

      hjetmass_qqbgg_bubble_pmmp_s124 = -ret/16q0*(0,1q0)
      return

      end function

      complex*32 function hjetmass_qqbgg_bubble_pmmp_s12
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt

      t1 = za(i2, i3)
      t2 = za(i1, i3)
      t3 = zb(i3, i1)
      t4 = zb(i3, i2)
      t5 = 0.1q1 / t4
      t6 = t2 * t3
      t7 = t6 * t5 + t1
      t8 = zb(i4, i1)
      t9 = za(i1, i4)
      t10 = za(i3, i4)
      t11 = 0.1q1 / t2
      t12 = t3 * t5
      t13 = t1 * t11
      t14 = t13 + t12
      t15 = za(i2, i4)
      t16 = 0.1q1 / t9
      t17 = t15 * t16
      t18 = t17 + t12
      t19 = zb(i4, i3)
      t20 = zb(i4, i2)
      t21 = 0.1q1 / t20
      t22 = t21 * t8 + t17
      t23 = t15 * t2
      t24 = t23 * t16
      t25 = t24 - t1
      t26 = t15 * t20
      t27 = t9 * t3
      t28 = t15 * t4
      t29 = t27 + t28
      t30 = t29 * t1
      t31 = t23 * (t28 * t16 + t3)
      t32 = t30 - t31
      t33 = t3 * t20
      t34 = t33 * t5 - t8
      t30 = t30 + t6 * (t27 * t5 + t15)
      t35 = t9 * t20
      t36 = t2 * t4
      t37 = t35 + t36
      t38 = t1 * t3
      t39 = t15 * t8
      t40 = t39 + t38
      t41 = t1 * t4
      t42 = t9 * t8
      t43 = -t41 + t42 - t26 + t6
      t43 = t43 ** 2
      t44 = 4 * t40 * t37 + t43
      t44 = sqrt(t44)
      t45 = -t41 + t42 - t26 + t6 + t44
      t46 = 0.1q1 / t37
      t47 = (0.1q1 / 0.2q1)
      t48 = t47 * t45 * t46
      t49 = t48 - t12
      t44 = t41 - t42 + t26 - t6 + t44
      t50 = t47 * t44 * t46
      t51 = -t50 - t12
      t52 = t1 ** 2
      t53 = t3 ** 2
      t54 = t41 * t26
      t55 = t52 * t4 ** 2
      t56 = t2 * t53
      t57 = t9 * t1
      t58 = -t23 + t57
      t59 = -t41 - t6
      t60 = t10 * t19
      t61 = t4 * t8
      t62 = 2
      t63 = -t48 * t20 + t8
      t64 = 0.1q1 / t37
      t65 = t35 * t1
      t66 = t59 * t2
      t67 = t6 * (t42 + t6)
      t68 = t2 * t45 * t64
      t69 = t9 * (t68 - t1) + t23
      t68 = t68 + t1
      t37 = 4 * t40 * t37 + t43
      t37 = sqrt(t37)
      t43 = -t41 + t42 - t26 + t6 + t37
      t43 = 0.1q1 / t43
      t70 = t62 * t40
      t71 = -t70 * t43 - t48
      t72 = t46 * (t45 + t44)
      t73 = -t48 * t2 - t1
      t74 = -t41 + t6
      t75 = t64 ** 2
      t76 = t46 * (t60 * (-t23 - t57) + t57 * (t41 - t6) + t23 * t74)
      t77 = t42 * t58
      t78 = t26 * t58
      t79 = t1 * (t15 * (t60 - t6) + t27 * t1)
      t80 = t9 * (-t41 + t60) + t23 * t4
      t81 = t75 * t45 ** 2 * t2
      t82 = (0.1q1 / 0.4q1)
      t83 = t62 * t39 * t58
      t84 = -t47 * t45 * (t35 * (-t23 + t57) * t75 * t45 + t76) + t82 * 
     #t81 * t80 + t79 + (t77 - t78) * t64 * t45 + t83
      t48 = -t48 - t17
      t85 = t2 ** 2
      t86 = t1 * (t15 * (t61 + t33) + t60 * t3)
      t87 = t52 * t4
      t88 = t61 + t33
      t46 = t46 * (-t23 * t88 + t57 * t88 + t60 * t74)
      t74 = t60 * t4 + t88 * t9
      t88 = t62 * t3 * (t2 * (t39 + t38) + t87)
      t53 = -t47 * t45 * (t4 * (t2 * (-t41 - t6) - t65) * t75 * t45 + t4
     #6) + t82 * t81 * t74 - t86 - (t53 * t85 - t41 * (t41 + t26) + t42 
     #* t6) * t64 * t45 - t88
      t81 = t16 ** 2
      t89 = t16 * t81
      t90 = t15 ** 2
      t3 = -t62 * t31 * t1 + t29 * t52 + t85 * (t4 * t15 * t90 * t81 + t
     #3 * t90 * t16)
      t29 = -t50 + t17
      t31 = t5 ** 2
      t91 = t5 * t31
      t92 = t30 * t8 + t33 * (-t1 * t15 + t12 * (-t23 - t57) - t56 * t9 
     #* t31)
      t93 = t50 * t20 + t8
      t94 = t2 * t44 * t64
      t95 = t9 * (t94 + t1) - t23
      t94 = t94 - t1
      t37 = t41 - t42 + t26 - t6 + t37
      t37 = 0.1q1 / t37
      t70 = t70 * t37 + t50
      t50 = t50 * t2 - t1
      t96 = t75 * t44 ** 2 * t2
      t76 = t47 * t44 * (t35 * (t23 - t57) * t75 * t44 + t76) + t82 * t9
     #6 * t80 + t83 + t79 + (-t77 + t78) * t64 * t44
      t46 = t47 * t44 * (t4 * (t2 * (t41 + t6) + t65) * t75 * t44 + t46)
     # + t82 * t96 * t74 - t88 - t86 - (t54 + t55 - t67) * t64 * t44
      t27 = t2 * (t27 - t28) + t41 * t62 * t9
      t28 = t60 + t6
      t42 = t42 + t26
      t47 = -t41 + t6
      t74 = t57 * t47
      t47 = t23 * t47
      t79 = t36 * t58 * t64
      t80 = t35 * t58 * t64
      t82 = 0.1q1 / t51
      t83 = 0.1q1 / t29
      t22 = 0.1q1 / t22
      t10 = 0.1q1 / t10
      t86 = 0.1q1 / t49
      t19 = 0.1q1 / t19
      t88 = 0.1q1 / t18
      t96 = -0.1q1 / t48
      t18 = 0.1q1 / t18
      t14 = 0.1q1 / t14
      t97 = t34 ** 2
      t98 = (t26 * t16 + t8) * t21
      t99 = t98 * t22
      t100 = t7 * t11
      t101 = t34 * t14
      t10 = t10 * t19
      t19 = 0.1q1 / t72
      t70 = 0.1q1 / t70
      t49 = -0.1q1 / t49
      t71 = 0.1q1 / t71
      t51 = -0.1q1 / t51
      t29 = -0.1q1 / t29
      t48 = 0.1q1 / t48
      t72 = t49 ** 2
      t102 = t71 ** 2
      t103 = t51 ** 2
      t104 = t29 ** 2
      t105 = t43 ** 2
      t106 = t63 ** 2
      t107 = t48 ** 2
      t108 = t63 * t73
      t109 = t53 * t72 * t31
      t110 = t73 * t84 * t81 * t107
      t111 = t50 * t93
      t112 = t93 ** 2 * t31 * t103
      t113 = t70 * t37
      t114 = t113 * t44
      t47 = (t114 * (t112 * (t46 * t70 + t62 * (t54 + t55 + t4 * (t66 - 
     #t65) * t64 * t44 - t67) - t61 * t95 - t60 * (t4 * t94 + t6) - t33 
     #* t95) + (t76 * (t2 * t93 + t50 * (-t70 * t93 + t20)) + t111 * (t6
     #2 * (-t80 * t44 - t77 + t78) + t47 - t79 * t44 - t60 * (-t9 * t94 
     #+ t23) - t74)) * t104 * t81) + (t71 * (-t106 * (-t62 * (t54 + t55 
     #+ t4 * (-t66 + t65) * t64 * t45 - t67) - t61 * t69 - t60 * (t4 * t
     #68 - t6) - t33 * t69) * t72 * t31 + (t84 * (t2 * t63 + t20 * t73) 
     #+ t108 * (t62 * (t80 * t45 - t77 + t78) - t60 * (t68 * t9 + t23) +
     # t47 + t79 * t45 - t74)) * t107 * t81) + t63 * (t109 * t63 - t110)
     # * t102) * t43 * t45) * t75 * t19
      t65 = t10 * t64
      t66 = t65 * t19 * t40
      t67 = t19 ** 2
      t68 = t64 * t45
      t44 = t44 * t64
      t43 = t71 * t43
      t69 = t43 * t63
      t71 = t10 * t75
      t44 = t71 * t67 * t40 * (t69 * (t53 * (t68 * t63 * t31 * t49 * t72
     # + (-t68 * t20 + t63) * t31 * t72) - t110 * (t68 * t48 + 1)) + t11
     #3 * t112 * t46 * (t44 * t51 - 1) + t113 * (-t44 * t20 * t46 * t31 
     #* t103 + (-t44 * t29 * t104 + t104) * t81 * t76 * t50) * t93)
      t53 = t18 ** 2
      t68 = t96 ** 2
      t72 = t83 ** 2
      t74 = t32 * t31
      t77 = t25 * t3
      t78 = t77 * t75 * t68 * t72
      t79 = t74 * t53 + t78
      t80 = t82 ** 2
      t84 = t86 ** 2
      t94 = t88 ** 2
      t95 = t92 * t14
      t107 = t97 * t80 * t84 * t75
      t112 = t30 * t81
      t3 = t10 * (t8 * (t98 * t15 * t89 * t79 * t22 ** 2 + t81 * (-t17 *
     # t79 + t98 * ((t17 * (-t2 * t3 - t25 * (t23 * t28 - t57 * t28 + t6
     #2 * (t85 * t90 * t4 * t16 + t23 * t42 - t57 * t42) + t87 * t9 - 3 
     #* t41 * t23)) - t77) * t68 * t75 * t72 - t31 * t53 * (t17 * t27 + 
     #t32))) * t22) + t13 * t14 * t31 * (-t34 * (t34 * t92 * t75 * t84 *
     # t80 + t112 * t94) + t12 * (t107 * (t95 + t62 * (t56 * (t35 * t5 +
     # t2) + t54 + t55) + t26 * t6 - t61 * t58 - t60 * t59 - t38 * (-4 *
     # t36 - 3 * t35)) + (-t27 * t34 + t30 * (t101 - t20)) * t94 * t81))
     #)
      t4 = t66 * t5 * (t43 * t1 * t106 * t49 + t113 * t93 * t51 * (t1 * 
     #t93 + t50 * t8) + t43 * t108 * t8 * t49)
      t6 = t114 * t93
      t9 = t10 * t64 * t75 * t19 * t67 * t40 * (t45 * (t43 * t109 * t106
     # - t69 * t110) + t6 * (-t93 * t46 * t31 * t103 + t76 * t50 * t81 *
     # t104))
      t2 = t1 * t20 + t2 * t8
      t15 = t10 * (t38 * (t101 * t91 * (-t20 * t92 * t80 * t84 * t75 + t
     #112 * t88 * t94) * t11 - t95 * t107 * t91 * (t82 + t86) * t11) + t
     #99 * t39 * t89 * (t74 * t18 * t53 + t78 * (t83 + t96)))
      t20 = t65 * t101 * t38 * t31 * t86 * t82 * (t13 * t20 + t8)
      ret = -8 * t10 * (-t52 * t7 * t97 * t11 ** 2 * t14 ** 2 * t5 * t64
     # * t86 * t82 + t101 * t52 * t11 * t16 * t88 * t5 * (t100 * t14 + 1
     #) + t5 * t16 * t1 * (t22 * t18 * (t98 + t17) + (t100 + t12) * t14 
     #* t88) * t8 + t16 * t21 * t22 * (t98 * t25 ** 2 * t83 * t22 * t64 
     #* t96 + (t25 * (-t99 + 1) + t24) * t18 * t5) * t8 ** 2 + t33 * t52
     # * t11 * t16 * t14 * t88 * t31) + 64 * t66 * (t40 * (t73 ** 2 * t6
     #3 * t16 * t105 * t48 * t102 - t111 * t37 ** 2 * t70 ** 2 * (t50 * 
     #t16 * t29 + t93 * t5 * t51) + t73 * t106 * t5 * t105 * t49 * t102)
     # + t47) + 128 * t44 - 16 * t3 + 12 * t13 * t65 * t101 * t5 * t86 *
     # t82 * (t1 * t34 + t7 * t8) + 48 * t4 - 256 * t9 - 40 * t71 * t5 *
     # t19 * t40 * (t69 * t2 * t49 * t45 - t6 * t51 * t2) - 32 * t15 + 2
     #0 * t20

      hjetmass_qqbgg_bubble_pmmp_s12 = -ret/16q0*(0,1q0)
      return

      end function

      complex*32 function hjetmass_qqbgg_bubble_pmmp_s34
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt

      t1 = za(i2, i3)
      t2 = zb(i3, i1)
      t3 = za(i1, i3)
      t4 = zb(i4, i1)
      t5 = zb(i4, i2)
      t6 = t3 * t4
      t7 = t1 * t5
      t8 = t7 + t6
      t9 = za(i1, i2)
      t10 = zb(i2, i1)
      t11 = zb(i3, i2)
      t12 = za(i1, i4)
      t13 = za(i2, i4)
      t14 = t12 * t2
      t15 = t13 * t11
      t16 = t14 + t15
      t17 = t3 * t2
      t18 = t1 * t11
      t19 = t12 * t4
      t5 = t13 * t5
      t20 = -t19 - t5 + t17 + t18
      t20 = t20 ** 2
      t21 = 4 * t16 * t8 + t20
      t21 = sqrt(t21)
      t22 = t21 - t19 - t5 + t17 + t18
      t21 = -t21 - t19 - t5 + t17 + t18
      t23 = t9 * t10
      t24 = (-t23 + t19 + t5) * t13
      t25 = t13 * (t18 + t17)
      t26 = t24 - 2 * t25
      t27 = 0.1q1 / t8
      t28 = t27 ** 2
      t29 = t27 * t28
      t30 = 2 * t4
      t31 = t13 * t16
      t32 = 0.1q1 / t22
      t33 = 0.1q1 / t21
      t34 = t32 ** 2
      t35 = t32 * t34
      t36 = t33 ** 2
      t37 = t33 * t36
      t38 = t8 ** 2
      t39 = t38 ** 2
      t40 = t38 * t39
      t41 = t31 * t4
      t42 = t8 * t39 * t34 * t36 * (t32 + t33)
      t31 = t2 * t31 * (384 * t40 * t34 * t36 * (t36 + t34) + 512 * t40 
     #* t35 * t37) * t28
      t39 = 16 * t39
      t40 = t23 + t19
      t43 = t1 * t40
      t44 = 2 * t13 * (t19 + t5) + (t23 - t17 - t18) * t13
      t45 = (t18 + t17) * t27
      t46 = 0.1q1 / t8
      t47 = (-t23 + t19 + t5) * t46
      t48 = (t7 + t6) * t28
      t49 = t21 * t13
      t50 = t23 + t17 + t18
      t51 = t21 ** 2
      t52 = t1 * t51 * t28
      t53 = (0.1q1 / 0.2q1)
      t54 = (0.1q1 / 0.4q1)
      t55 = -t1 * t21 * t51 * t29 * t8 / 8
      t56 = -t53 * t49 * (t48 * t21 + t47) + t54 * t52 * t50 + t13 * (t4
     #5 * t21 + t14 + t15) + t55
      t20 = 4 * t16 * t8 + t20
      t20 = sqrt(t20)
      t57 = -t20 - t19 - t5 + t17 + t18
      t57 = 0.1q1 / t57
      t58 = 2 * t16
      t59 = t53 * t21 * t46
      t60 = -t58 * t57 - t59
      t61 = t46 * (-t22 + t21)
      t62 = -t59 * t4 + t2
      t63 = -t53 * t4 * t22 * t46 + t2
      t64 = (-t5 - t19) * t27
      t65 = t22 ** 2
      t66 = t22 * t13
      t29 = -t1 * t22 * t65 * t29 * t8 / 8
      t67 = (0.3q1 / 0.4q1) * t7 * t13
      t68 = -t53 * t66 * (t46 * (t23 - t17 - t18) + t6 * t22 * t28) - t5
     #4 * t43 * t65 * t28 + t13 * (t64 * t22 + t14 + t15) + t29 - t67 * 
     #t65 * t28
      t20 = t20 - t19 - t5 + t17 + t18
      t20 = 0.1q1 / t20
      t69 = t53 * t22 * t46
      t58 = -t58 * t20 - t69
      t65 = t1 * t65 * t28
      t29 = -t53 * t66 * (t48 * t22 + t47) + t54 * t65 * t50 + t29 + t13
     # * (t45 * t22 + t14 + t15)
      t45 = t1 * t22 * t27
      t47 = (0.3q1 / 0.4q1) * t65 * t8
      t48 = -t23 - t17
      t11 = t1 ** 2 * t11
      t50 = (-t7 - t6) * t27
      t65 = t1 * t21 * t27
      t70 = (0.3q1 / 0.4q1) * t52 * t8
      t14 = t53 * t49 * (t46 * (-t23 + t17 + t18) - t6 * t21 * t28) - t5
     #4 * t52 * t40 + t55 + t13 * (t64 * t21 + t14 + t15) - t67 * t51 * 
     #t28
      t9 = 0.1q1 / t9
      t10 = 0.1q1 / t10
      t15 = t32 + t33
      t40 = t13 * t4
      t46 = t1 * t2
      t51 = 0.1q1 / t61
      t52 = 0.1q1 / t60
      t53 = 0.1q1 / t58
      t54 = t46 + t40
      t55 = t57 * t36 * t52
      t58 = t9 * t10
      t60 = t58 * t51 * t16 * t8
      t61 = -t56 + t14
      t64 = t51 ** 2
      t47 = t63 * (-t68 + t29) * t20
      t14 = (t56 - t14) * t52
      t37 = t37 * t52 * t57
      t52 = t58 * t64 * t16
      t3 = t52 * t8 * (t35 * (t47 * t53 ** 2 + (t4 * (t68 - t29) + t63 *
     # (2 * t13 * (t4 * (t3 * t22 * t27 + t12) + t5) + t19 * t45 - t25 +
     # t23 * (t45 + t13) + 3 * t66 * t7 * t27 + t11 * t22 * t27 - t45 * 
     #t48 - t24 + 2 * t13 * (t50 * t22 + t17 + t18))) * t20 * t53) + t37
     # * (t4 * t61 + t62 * (t14 + 2 * t13 * (t4 * (t3 * t21 * t27 + t12)
     # + t5) + t65 * t19 - t25 + t23 * (t65 + t13) + 3 * t49 * t7 * t27 
     #+ t11 * t21 * t27 - t65 * t48 - t24 + 2 * t13 * (t50 * t21 + t17 +
     # t18))))
      ret = -64 * t60 * (t54 * t20 * t53 * t34 - t55 * t54) + 192 * t60 
     #* (t55 * (t1 * t62 - t4 * (-t59 * t1 - t13)) + t20 * t34 * t53 * (
     #-t1 * t63 + t4 * (-t69 * t1 - t13))) + 1024 * t3 + 4096 * t58 * t5
     #1 * t64 * t16 * t8 * (t47 * t35 * t53 + t37 * t61 * t62) - 6144 * 
     #t52 * t38 * (t14 * t36 ** 2 * t57 * t62 + t47 * t34 ** 2 * t53) - 
     #8 * t58 * (-128 * t28 * (t2 * t26 + t41) * t42 + t39 * t28 * (t2 *
     # (2 * (t23 + t17 + t18) * t1 - 4 * t13 * t8) + t30 * t26) * t34 * 
     #t36 + 128 * t28 * (t2 * t44 + t41) * t42 - t39 * t28 * (t2 * (t13 
     #* (-6 * t7 - 4 * t6) - 2 * t43) + t30 * t44) * t34 * t36) - 128 * 
     #t33 * t10 * t9 * t32 * t38 * (t40 * t15 + t46 * t15)

      hjetmass_qqbgg_bubble_pmmp_s34 = -ret/16q0*(0,1q0)
      return

      end function

      complex*32 function hjetmass_qqbgg_bubble_pmpm_mhsq
     &     (i1,i2,i3,i4,za,zb,mt,p,flip)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          double precision mt
          real*16 p(mxpart,4)
          complex*32 hjetmass_qqbgg_bubble_pmpm_s124
          complex*32 hjetmass_qqbgg_bubble_pmpm_s123
          complex*32 hjetmass_qqbgg_bubble_pmpm_s34
          complex*32 hjetmass_qqbgg_bubble_pmpm_s12
          logical flip



      hjetmass_qqbgg_bubble_pmpm_mhsq = 
     & -hjetmass_qqbgg_bubble_pmpm_s124(i1,i2,i3,i4,za,zb,mt,p,flip)
     & -hjetmass_qqbgg_bubble_pmpm_s123(i1,i2,i3,i4,za,zb,mt,p,flip)
     & -hjetmass_qqbgg_bubble_pmpm_s34(i1,i2,i3,i4,za,zb,mt)
     & -hjetmass_qqbgg_bubble_pmpm_s12(i1,i2,i3,i4,za,zb,mt)
      return

      end function

      complex*32 function hjetmass_qqbgg_bubble_pmpm_s123
     &     (i1,i2,i3,i4,za,zb,mt,p,flip)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt
          real*16 p(mxpart,4)
          complex*32 alpha
          logical flip

      alpha = (za(i1,i2)*zb(i2,i1) + za(i1,i3)*zb(i3,i1) +
     & za(i2,i3)*zb(i3,i2))/(za(i1,i2)*zb(i2,i1) + za(i1,i3)*zb(i3,i1))

      p(5,:) = real(alpha)*p(i1,:)
      p(6,:) = (1q0-real(alpha))*p(i1,:) + p(i2,:) + p(i3,:)
      if (flip .eqv. .true.) then
          call spinoru_qp(6,p,zb,za)
      else
          call spinoru_qp(6,p,za,zb)
      end if

      t1 = za(6, i4)
      t2 = zb(i3, 5)
      t3 = zb(i3, i1)
      t4 = za(5, i2)
      t5 = zb(i2, 6)
      t6 = za(5, i3)
      t7 = zb(i3, 6)
      t8 = za(5, i4)
      t9 = zb(i4, 6)
      t10 = t8 * t9
      t11 = t4 * t5
      t12 = t6 * t7
      t13 = t12 + t10 + t11
      t14 = zb(i2, i1)
      t15 = za(i1, i2)
      t16 = za(i2, i3)
      t17 = zb(i3, i2)
      t18 = t15 * t14
      t19 = za(i1, i3) * t3
      t20 = t16 * t17
      t21 = t20 + t18 + t19
      t22 = za(6, i1)
      t23 = zb(i1, 6)
      t24 = zb(i2, 5)
      t25 = za(6, i2)
      t26 = za(6, i3)
      t27 = zb(i4, 5)
      t28 = t25 * t24
      t29 = t26 * t2
      t30 = t1 * t27
      t31 = t30 + t28 + t29
      t32 = t22 * t23
      t33 = t4 * t24
      t34 = t25 * t5
      t35 = t6 * t2
      t36 = t26 * t7
      t37 = t8 * t27
      t38 = t1 * t9
      t39 = -t37 + t38 + t36 + t34 - t35 + t32 - t33
      t39 = t39 ** 2
      t40 = 4 * t31 * t13 + t39
      t40 = sqrt(t40)
      t41 = -t37 + t38 + t36 + t34 - t35 + t32 - t33 - t40
      t40 = -t37 + t38 + t36 + t34 - t35 + t32 - t33 + t40
      t42 = zb(i4, i3)
      t43 = t30 + t28 + t29
      t44 = t1 * t25
      t45 = t44 * t43
      t46 = za(i2, i4)
      t47 = za(i1, i4)
      t48 = za(i3, i4)
      t49 = t47 * zb(i4, i1)
      t50 = t42 * t48
      t51 = t46 * zb(i4, i2)
      t52 = t51 + t18 + t19 + t20 + t32 + t34 + t36 + t38 + t49 + t50
      t53 = t34 + t32 - t33
      t54 = t34 * t33
      t53 = t53 ** 2 + 4 * t54
      t53 = sqrt(t53)
      t55 = t34 + t32 - t33 + t53
      t56 = 0.1q1 / t5
      t57 = (0.1q1 / 0.2q1)
      t58 = t55 * t56
      t59 = t58 * t57
      t60 = t59 - t25
      t61 = 0.1q1 / t4
      t62 = t59 * t7 * t61 + t2
      t63 = 0.1q1 / t8
      t64 = t1 * t63
      t59 = t59 * t61
      t65 = t59 - t64
      t53 = t34 + t32 - t33 - t53
      t66 = t61 * t56
      t67 = t66 * (t55 - t53)
      t68 = 0.1q1 / t55
      t69 = 2 * t28
      t70 = t69 * t68 + t59
      t71 = 0.1q1 / t9
      t72 = t27 * t71
      t73 = t59 + t72
      t74 = t53 * t56
      t75 = t74 * t57 - t25
      t76 = t66 * t57 * t53
      t77 = t76 * t7 + t2
      t78 = 0.1q1 / t53
      t69 = t69 * t78 + t76
      t79 = 0.1q1 / t13
      t80 = t41 * t79
      t81 = t66 * t53
      t82 = -t81 + t80
      t83 = t40 * t79
      t84 = -t81 + t83
      t85 = t76 - t64
      t86 = t76 + t72
      t87 = t58 * t61
      t88 = -t87 + t80
      t89 = -t87 + t83
      t83 = t83 * t57
      t90 = t83 * t4 - t25
      t91 = t83 * t8 - t1
      t92 = t83 * t7 + t2
      t39 = 4 * t31 * t13 + t39
      t39 = sqrt(t39)
      t93 = -t37 + t38 + t36 + t34 - t35 + t32 - t33 + t39
      t93 = 0.1q1 / t93
      t94 = 2 * t31
      t95 = t94 * t93 + t83
      t96 = t79 * (t41 - t40)
      t80 = t80 * t57
      t97 = t80 * t7 + t2
      t98 = 0.1q1 / t13
      t99 = t8 * t25
      t100 = t1 * t4
      t101 = t44 * (t51 + t50 + t49 + t36 + t32 + t20 + t18 + t19)
      t102 = -t99 - t100
      t103 = t38 + t34
      t104 = t8 ** 2
      t105 = t20 * t25
      t106 = t2 * t4
      t107 = t106 * t26
      t108 = t104 * t25
      t109 = t108 * t27
      t110 = t44 * t103
      t111 = t12 + t10 + t11
      t112 = t111 * t98
      t113 = t12 * t25
      t114 = t12 + t11
      t115 = t98 ** 2
      t116 = t115 ** 2
      t117 = t98 * t116
      t118 = t98 * t115
      t119 = t41 ** 2
      t120 = t119 ** 2
      t121 = t41 * t119
      t122 = t100 * t114
      t123 = t108 * t9
      t124 = (t4 * (t51 + t18 + t19 + t20 - t32 - t36 + t49 + t50) + t11
     #3) * t8
      t125 = t115 * t119
      t126 = (0.3q1 / 0.4q1)
      t127 = t126 * t125 * (t124 + t122 + t123) + (t100 * (-t51 - t50 - 
     #t49 - t20 - t18) + t99 * (-t51 - t50 - t49 + t35 - t18)) * t98 * t
     #41 + t101 + (t8 * (-t107 - t105) + t109 + t19 * t102 + t100 * (t35
     # + t33)) * t98 * t41 + t110 - 2 * t44 * (t112 * t41 + t33 + t35 + 
     #t37)
      t39 = -t37 + t38 + t36 + t34 - t35 + t32 - t33 - t39
      t39 = 0.1q1 / t39
      t94 = t94 * t39 + t80
      t128 = (-t37 - t33 - t35) * t98
      t129 = -t12 - t10 - t11
      t130 = t52 * t79
      t131 = t44 * t41
      t124 = t122 * t118 + t123 * t118 + t124 * t118
      t132 = t4 * t26
      t133 = t6 * t25
      t134 = -t133 + t132
      t135 = t99 + t100
      t136 = t35 + t33
      t137 = t18 * t135
      t138 = t2 * t8
      t139 = -t100 * t136 + t138 * t134 + t19 * t135 + t20 * t135 + t49 
     #* t135 + t50 * t135 + t51 * t135 - t109 + t137
      t140 = (0.1q1 / 0.4q1)
      t141 = (0.1q1 / 0.8q1)
      t142 = t57 * t131 * (t129 * t115 * t41 + t130) + t141 * t121 * t12
     #4 - t140 * t125 * t139 + t44 * (t128 * t41 + t28 + t29 + t30)
      t143 = -t51 - t19 - t20 - t49 - t50
      t144 = t37 + t35
      t145 = t4 ** 2
      t146 = t138 * t26
      t147 = t145 * t24
      t148 = t147 * t1
      t149 = t40 ** 2
      t150 = t149 ** 2
      t151 = t40 * t149
      t152 = t115 * t149
      t101 = -t126 * t152 * (t8 * (t4 * (-t51 - t18 - t19 - t20 + t32 + 
     #t36 - t49 - t50) - t113) - t122 - t123) + (t99 * t143 + t100 * (-t
     #51 - t49 + t35 - t20 - t19)) * t98 * t40 + t101 + (t4 * (-t50 * t1
     # - t146) + t148 + t18 * t102 + t99 * t144) * t98 * t40 + t110 - 2 
     #* t44 * (t112 * t40 + t33 + t35 + t37)
      t110 = t44 * t40
      t111 = -t140 * t152 * t139 + t141 * t151 * t124 - t57 * t110 * (t1
     #11 * t115 * t40 - t130) + t44 * (t128 * t40 + t28 + t29 + t30)
      t112 = t80 * t4 - t25
      t124 = t80 * t8 - t1
      t128 = t100 * t63 - t25
      t130 = t1 * t7
      t139 = t130 * t63 + t2
      t153 = t72 + t64
      t154 = t25 * t61
      t155 = t56 * (-t61 * (t32 + t18) + t24) - t154
      t156 = t25 * t53
      t157 = t56 * (-t32 + t18) + t25
      t158 = t53 ** 2
      t159 = t158 ** 2
      t160 = t53 * t158
      t161 = t66 * t158
      t162 = t25 ** 2
      t163 = t162 * t24
      t164 = t140 * t161 * t157 + t57 * t156 * t155 - t163
      t165 = t74 + t25
      t166 = t25 * t55
      t167 = t55 ** 2
      t168 = t167 ** 2
      t169 = t55 * t167
      t170 = t66 * t167
      t155 = t140 * t170 * t157 + t57 * t166 * t155 - t163
      t157 = t58 + t25
      t171 = -t37 + t18 + t32 + t34 - t35
      t172 = t43 * t25
      t173 = t171 * t25 * t98
      t174 = t18 * t125
      t175 = t174 * t4
      t176 = t44 * t24
      t177 = t12 * t1
      t178 = t8 * (-t36 + t35) - t177
      t179 = t5 * t8
      t180 = t179 * t162
      t181 = -t180 + t148
      t182 = t6 * t1
      t183 = t8 * t26
      t184 = t183 + t182
      t185 = t35 * t25
      t186 = t7 * t25
      t187 = t4 * t41 * t98
      t188 = t36 * t34
      t189 = t106 * t184
      t190 = t51 * t181
      t191 = t118 * t121
      t192 = t24 ** 2
      t193 = t1 ** 2
      t194 = t100 * t5
      t195 = t20 * t181
      t196 = t49 * t181
      t197 = t4 * t145 * t1
      t198 = t197 * t192
      t199 = t193 * t9
      t200 = t199 * t34 * t4
      t201 = t5 * t27
      t202 = t9 * t24
      t203 = t202 + t201
      t204 = t8 * t41 * t98
      t205 = t104 * t162
      t206 = t205 * t203
      t207 = t23 ** 2 * t22 ** 2
      t208 = t207 * t4
      t197 = t197 * t24 * t5
      t134 = t191 * (-t197 * t41 * t98 + t206 + t32 * (-t99 * (t51 + t50
     #) + t100 * (-t204 * t9 - t50 - t51)) + t187 * t34 * t8 * (t51 + t3
     #8) + t208 * (-t204 + t1)) + t191 * (-t180 * t50 + t195 + t196 - t1
     #98 - t200 + t18 * (-t100 * (t38 + t36) + t109) + t32 * (t100 * (-t
     #49 + t38) + t109 + t99 * (-t49 - t19)) - t109 * t33 + t187 * (t8 *
     # (t34 * (t50 + t49) + t32 * (-t36 + t33 - t34)) + t18 * (t8 * (t38
     # + t32 - t33) + t194))) + t191 * (t8 * (t24 * (t106 * t134 + t186 
     #* (t133 - t132)) + t34 * (t4 * (t20 * t41 * t98 - t29) + t185)) + 
     #t190 + t100 * (t33 * (t50 - t35) - t188) + t19 * (t179 * (t187 * t
     #25 - t162) + t148) + t32 * (t99 * (t35 - t20) + t100 * (t36 - t19 
     #- t20)) + t18 * (t178 * t25 + t189))
      t204 = t4 * t27
      t209 = t25 * t9
      t210 = -t209 + t204
      t211 = t8 * t24
      t212 = t1 * t5
      t213 = -t212 + t211
      t214 = -t138 - t130
      t215 = -t138 + t130
      t216 = t37 * t100
      t217 = t20 * t213
      t218 = t8 * t203
      t219 = t19 * t213
      t220 = t8 * t210
      t221 = t220 * t18 * t1
      t222 = t4 * t25
      t223 = t32 * t1
      t224 = t100 * t2
      t225 = t35 * t44
      t226 = t209 - t204
      t227 = -t99 - t100
      t228 = t212 - t211
      t229 = t28 * t8
      t230 = t228 * t25
      t231 = t230 * t50
      t232 = t18 * t227
      t233 = t107 * t8
      t234 = t12 * t44
      t235 = t19 * t25
      t236 = t18 * t193
      t237 = t51 * t162
      t238 = t51 + t49 + t50
      t239 = t49 + t20
      t240 = t5 ** 2
      t241 = t162 * t240
      t242 = t99 * t9
      t243 = t33 * t8
      t244 = t100 * t27
      t245 = t207 * t44
      t238 = t125 * (t25 * (t1 * (t34 * t239 + t241) + t229 * (-t49 + t3
     #5 + t37) + t226 * t5 * t193) + t245 + t223 * (t238 * t25 + t244) +
     # t175 * t34 * t8 + t25 * (t243 * t238 + t212 * (-t238 * t4 - t242)
     #) * t98 * t41) + t125 * (t25 * (t1 * (t38 * (-t33 + t32) + t188) +
     # t231 + t20 * (-t229 + t223)) + (-t177 * t5 * t162 + t32 * (-t234 
     #+ t232 - t233) + t18 * t181) * t98 * t41 + t235 * (t223 + t230) + 
     #t236 * t226 + t237 * t228) + t125 * (t26 * (t222 * (-t212 * t2 + t
     #214 * t24) + t223 * (t186 + t106)) + t18 * (t26 * (t215 * t25 - t2
     #24) - t225) + (t32 * (t25 * (t8 * (-t38 - t33) - t194) - t216) + t
     #221 + t222 * (t218 * t1 + t217 + t219)) * t98 * t41)
      t246 = t6 * t24
      t247 = t5 * t26
      t248 = t12 + t10
      t249 = t211 * t51
      t250 = (t247 + t246) * t7
      t251 = t99 * t248
      t252 = t32 * t8
      t253 = t116 * t120
      t254 = t2 * t5
      t255 = t24 * t7
      t256 = t255 + t254
      t257 = -t255 - t254
      t258 = t187 * t5
      t259 = t257 * t26
      t260 = t49 * t24
      t261 = t256 * t6
      t262 = t50 * t145 * t213
      t263 = t18 * t4
      t264 = t263 * t8
      t265 = t12 + t11
      t266 = -t12 - t10
      t267 = t30 * (-t33 + t18 + t32)
      t268 = t79 * t1
      t269 = t263 * t213 * t118
      t270 = t25 * t41
      t271 = -t51 - t49 - t50
      t272 = t230 * t18
      t119 = t270 * (t44 * t79 * (t24 * (t38 + t19 + t36) + t212 * t27) 
     #+ (t1 * (t33 * t271 + t34 * (-t35 + t32)) + t272) * t115 * t41) + 
     #t270 * (t1 * (t25 * (-t179 * t27 + t24 * t266) + t32 * (-t37 + t18
     # - t33 - t35) + t33 * (t37 - t19 - t20) - t37 * t18) * t115 * t41 
     #+ t268 * (t28 * (t51 + t20 + t49 + t50) + t267 + t29 * (t34 + t18 
     #+ t32 - t33)) + t269 * t119)
      t273 = -t240 * t25 + t261
      t274 = t5 * t41
      t275 = t51 + t19 + t20 - t36 + t49 + t50
      t276 = t18 * t8
      t277 = t276 * t248
      t278 = t252 * t248
      t279 = t5 * (t8 * (t275 * t4 + t113) + t122 + t123) + t277 - t278
      t280 = t24 * t79
      t281 = t162 * t1
      t282 = t281 * t33
      t283 = (0.3q1 / 0.2q1)
      t284 = (0.3q1 / 0.8q1)
      t285 = -(0.3q1 / 0.16q2)
      t286 = (0.1q1 / 0.16q2)
      t287 = (0.1q1 / 0.32q2)
      t119 = -t126 * t125 * t33 * t44 * (t4 * (t274 * t98 + t24) + t35) 
     #- t140 * t238 - t141 * t134 - t57 * t119 + t286 * (t253 * (t8 * (t
     #145 * (t260 + t259) + t6 * (-t5 * t7 * t162 + t106 * t34)) + t262 
     #+ t264 * (t258 + t35) - t205 * t5 * t9 + t145 * (-t51 * t5 + t261)
     # * t1 + t252 * (t4 * (-t51 - t35 - t258 - t49 - t50) - t113)) + t2
     #53 * (t4 * (t4 * (-t212 * t49 + t249) + t217 * t4 + t219 * t4 + t9
     #9 * (t218 + t250)) + t18 * (t4 * (t8 * (t37 - t36) - t177) - t251)
     # + t252 * (t4 * (-t37 - t19 - t20) - t242))) + t285 * t11 * t253 *
     # t44 * t265 - t284 * t191 * t44 * t4 * t273 + t287 * t117 * t41 * 
     #t120 * t4 * t279 + t283 * t282 * t41 * (t274 * t115 + t280) - t176
     # * (t173 * t41 + t172 - t175)
      t120 = t33 - t32
      t134 = t38 + t18
      t175 = t247 * t130
      t238 = t246 * t138
      t258 = t50 * t25
      t274 = t105 * t213
      t288 = t49 * t25
      t289 = t193 * t27
      t290 = -t30 - t29
      t291 = t37 * t4
      t292 = t222 * t49
      t293 = t222 * t51
      t294 = t222 * t50
      t295 = t222 * t19
      t296 = t199 * t5
      t297 = t104 * t24 * t27
      t298 = (t138 - t130) * t25
      t137 = t152 * (t25 * (t25 * (-t297 - t296) - t218 * t100 * t98 * t
     #40 + t130 * t33 * t26) - t245 + t18 * (t26 * (t298 + t224) + t225 
     #+ t1 * t8 * t226 * t98 * t40) + t32 * (t44 * (-t36 - t19 - t20) + 
     #(t25 * (t8 * (t38 + t33) + t194) + t137 + t233) * t98 * t40)) + t1
     #52 * (t25 * (-t281 * t240 + t107 * (t212 + t211)) + (t295 * t228 +
     # t223 * (t291 + t113) + t292 * t228 + t293 * t228 + t294 * t228 + 
     #t105 * t4 * t228) * t98 * t40 + t223 * (t25 * t271 + t290 * t4) - 
     #t264 * t152 * t34) + t152 * (t25 * (t25 * (-t238 - t175) + t258 * 
     #t213 + t274 + t288 * t213 + t51 * t25 * t213 + t219 * t25 + t289 *
     # t11 + t199 * t120) + (-t148 * t18 + (t134 * t8 + t177) * t162 * t
     #5) * t98 * t40 + t236 * t210)
      t194 = t199 * t4
      t224 = t145 * t228
      t264 = t35 * t99
      t299 = t118 * t151
      t300 = -t51 - t19 - t20 - t49 - t50
      t301 = t4 * t8
      t302 = t38 * t18
      t303 = t35 * t24
      t181 = t299 * (t189 * t18 + t206 + t99 * t4 * (-t303 + t259) + t20
     #8 * (-t8 * t40 * t98 + t1) - t197 * t40 * t98 + t32 * (t36 * t100 
     #+ t264)) + t299 * (t18 * (-t36 * t135 + t298 * t6) + t32 * (t301 *
     # (-t38 + t18 - t34 - t36) * t98 * t40 + t100 * (-t51 - t19 - t20 -
     # t49) + t99 * t300) + t301 * (t34 * (t51 + t19 + t20 + t38 + t49 +
     # t50) + t302) * t98 * t40 - t11 * t36 * t44) + t299 * (t24 * (t12 
     #* t8 * t162 - t109 * t4 + t145 * (t183 - t182) * t2) + t190 + t195
     # + t196 - t198 + t18 * (t224 * t98 * t40 + t109 - t194) + t32 * (t
     #147 * t8 * t40 * t98 - t50 * t100 + t109 + t194) + t50 * t181 + t1
     #9 * t181 + t34 * (-t194 + t264))
      t189 = t209 + t204
      t190 = t12 * t162
      t195 = t32 * t104 * t189
      t196 = t219 * t145
      t197 = t217 * t145
      t116 = t116 * t150
      t198 = -t34 + t33
      t208 = t179 * t145
      t264 = t35 * t5
      t298 = t51 + t19 + t20 + t36 + t49 + t50
      t304 = t37 + t35
      t305 = t266 * t25
      t306 = t34 * t304
      t79 = t110 * (t79 * (t44 * t203 + t29 * (-t33 + t18)) + (t18 * (-t
     #37 + t32) + t33 * (t37 - t19 - t20 - t32) - t37 * t32) * t115 * t4
     #0) + t40 * t25 * (t269 * t149 + t268 * (t28 * t298 + t267 + t29 * 
     #(t34 + t32)) + (t1 * (t24 * (t271 * t4 + t305) - t306 + t32 * (-t3
     #5 + t34)) + t272) * t115 * t40)
      t149 = t5 * t40
      t79 = -t126 * t152 * t33 * t44 * (t4 * (t149 * t98 + t24) + t35) +
     # t140 * t137 - t141 * t181 + t283 * t282 * t40 * (t115 * t149 + t2
     #80) - t284 * t299 * t44 * t4 * t273 + t285 * t116 * t11 * t44 * t2
     #65 + t286 * (t116 * (t252 * (t4 * (-t11 * t40 * t98 - t19 - t20 - 
     #t35 - t49 - t50 - t51) - t113) + t99 * (t4 * (t264 + t250) + t10 *
     # t198) + t18 * (t208 * t40 * t98 + t178 * t4 - t251)) + t116 * (t8
     # * (t259 * t145 - t190 * t5) + t262 - t195 + t196 + t197 + t204 * 
     #(t34 + t18) * t104 + t261 * t145 * t1 + t49 * t145 * t213 + t51 * 
     #t145 * t213)) + t287 * t117 * t40 * t150 * t4 * t279 - t57 * t79 -
     # t176 * (-t263 * t152 + t173 * t40 + t172)
      t117 = t204 * t71 + t25
      t137 = t72 * t7 - t2
      t149 = t53 * t61
      t150 = t55 * t61
      t173 = t281 * t24
      t80 = t80 * t5 + t24
      t83 = t83 * t5 + t24
      t178 = -t58 + t25
      t181 = t61 ** 2
      t203 = t61 * t181
      t250 = t56 ** 2
      t259 = t56 * t250
      t262 = t297 * t56
      t265 = t38 * t8
      t267 = t10 * t87
      t268 = t1 * t178
      t269 = t18 * t250
      t271 = t269 * t61
      t272 = t32 * t250
      t280 = t61 * t167
      t282 = t58 - t25
      t307 = t186 * t61
      t308 = t51 * t61
      t309 = t19 * t61
      t310 = t211 * t56
      t311 = t25 * t56
      t312 = t58 * t24
      t313 = t18 * t56
      t314 = t209 * t61
      t315 = t314 - t27
      t316 = t26 * (t138 + t130)
      t317 = t193 * t315
      t318 = t18 * t1
      t319 = t44 * t9
      t320 = t313 * t1
      t321 = t32 * t56
      t322 = t170 * (t25 * (t24 * (t56 * (-t316 + t99 * (-t51 - t19 - t2
     #0 + t35) * t61) + t150 * t8 * (t51 + t19 + t20 + t38) * t250) + t3
     #17) + t320 * (t61 * (t37 * t58 + t319) - t29) + t321 * (t1 * (t319
     # * t61 + t29) + t87 * (t8 * (-t30 - t28 - t29) - t318)))
      t323 = t51 + t19 + t20 - t35 + t49 + t50
      t324 = (t36 - t35) * t8
      t325 = t28 * t27
      t326 = t99 * t323
      t327 = t326 * t61
      t328 = t154 * (t177 + t324)
      t329 = t24 * t56
      t330 = t329 * t1
      t331 = t250 * t169
      t332 = t87 * t8 - t1
      t333 = t38 * t61
      t334 = -t333 + t24
      t335 = t99 * t61
      t336 = t335 * t27
      t337 = t276 * t66
      t338 = t32 * t66
      t339 = t331 * t61
      t340 = t154 * t6
      t341 = (-t340 + t26) * t7
      t342 = t341 + t35
      t343 = t32 * t87
      t344 = t331 * t181
      t239 = t344 * (t25 * (t8 * (t329 * t342 - t336) + t199) + t207 * t
     #56 * t332 + t313 * (t1 * (t36 - t55 - t35) + t8 * (-t343 - t29))) 
     #+ t339 * (t25 * (-t51 * t58 * t8 * t181 + t316 * t61) + t329 * (t1
     # * (-t50 - t19 + t35) - t146) + t337 * (t58 * t334 - t336) + t338 
     #* (t1 * (t51 + t19 + t20 + t267 + t49 + t50) + t335 * (-t37 + t55)
     # + t36 * t332 - t211 * t58) + (t51 + t19 - t35 + t50) * t8 * t181 
     #* t162) + t331 * (t104 * (-t202 * t203 * t56 * t162 + t325 * t181 
     #* t56) + t8 * (t203 * t239 * t162 + t58 * t203 * (-t50 - t19 - t20
     # - t38 - t49) * t25) + t321 * t181 * (-t199 + t327) + t313 * t181 
     #* (t199 + t328) + t330 * (t61 * (-t51 + t55 - t20 - t49) + t24))
      t332 = t181 * (t251 * t56 + t8 * (t18 * t248 - t32 * t248) * t250)
     # + t1 + t66 * (t275 * t8 + t177)
      t336 = t259 * t181
      t345 = t56 * (-t37 + t18 + t32 - t35) + t25
      t346 = t176 * (t166 * t345 * t61 - t269 * t280 + t172)
      t347 = -t51 - t19 - t20 + t36 - t49 - t50
      t348 = t113 * t61
      t251 = t251 * t61
      t340 = t335 * (t7 * (t340 - t26) - t35)
      t349 = t252 * t66
      t350 = t336 * t168
      t242 = t242 * t61
      t351 = t109 * t61
      t352 = t323 * t1
      t353 = t18 * t66
      t205 = t350 * (t353 * (t8 * (t242 - t55 - t35 + t36) + t177) - t35
     #1 + t146 + t352 + t349 * (t242 + t55 + t19 + t20 + t35) + t205 * t
     #9 * t181) + t350 * (t337 * (t348 - t37) + t340 + t329 * (t347 * t8
     # - t177 - t251) + t349 * (t51 + t37 + t348 + t49 + t50))
      t242 = -t37 + t18 - t35
      t337 = t202 * t56 + t27
      t354 = t337 * t8
      t355 = t44 * t61
      t260 = t166 * (t24 * (t150 * (-t1 * (t32 + t19 + t20) + t154 * (-t
     #276 - t177)) * t250 + t1 * (t154 * (t36 + t19 + t20) - t30) * t56 
     #+ t336 * t276 * t167) + t355 * t29) + t44 * t55 * (t50 * t66 * t24
     # * t178 + t51 * t66 * t24 * t178 + t260 * t66 * t178 + t355 * t337
     # + t353 * (t150 * (-t37 * t56 + t25) + t29 + t30 - t170) + t338 * 
     #(t150 * (t242 * t56 + t25) + t29 + t30) + t87 * (t329 * t37 + t154
     # * (-t354 - t35)) - t329 * t29)
      t337 = t150 + t24
      t356 = t176 * t250
      t357 = t356 * t167 * (t61 * (t35 + t55) + t24)
      t358 = (t329 * t7 + t2) * t6
      t359 = -t358 * t56 + t25
      t360 = t44 * t56 * t181
      t361 = t12 * t66 + 1
      t205 = -t126 * t357 - t140 * (t322 + t170 * (t25 * (t1 * (t309 * t
     #178 - t29 + t66 * (-t32 * t55 + t37 * t55 + t207) + t307 * (-t87 *
     # t6 + t26) + t308 * t178) + t49 * t66 * (-t1 * t55 + t211 * t282) 
     #+ t50 * t61 * (t310 * t282 + t268)) + t313 * (t1 * (t312 - t30) + 
     #t8 * (-t55 * t181 * t162 + t311 * t167 * t181))) + t280 * (t25 * (
     #t25 * (-t265 * t55 * t250 * t181 + t66 * (t49 * t1 + t262)) + t281
     # * t61 + t268 * t20 * t66 + t271 * (t1 * (-t267 - t35 + t36) - t14
     #6) - t199 * t24 * t250) + t272 * (t25 * (t298 * t1 * t61 + t58 * (
     #t8 * (-t38 - t18) - t177) * t181) + t289))) + t141 * t239 + t283 *
     # t173 * t58 * t337 + t284 * t360 * t169 * t359 + t285 * t44 * t168
     # * t181 * t250 * t361 - t286 * t205 + t287 * t336 * t55 * t168 * t
     #332 - t57 * t260 - t346
      t239 = t71 ** 2
      t76 = t76 * t8 - t1
      t59 = t59 * t8 - t1
      t260 = t61 * t158
      t267 = t74 - t25
      t268 = -t74 + t25
      t282 = t53 * t250
      t322 = t199 * t25
      t346 = t1 * t268
      t357 = t74 * t24
      t362 = t156 * t61
      t307 = t161 * (t25 * (t346 * t309 + t317 + t66 * (t325 * t104 + (-
     #t49 + t37 - t348) * t53 * t1) + t1 * (t307 - t2) * t26 + t346 * t3
     #08 + t346 * t20 * t61) + t223 * t56 * (-t362 + t30) + t245 * t66 +
     # t313 * (t1 * (t357 - t30) + t8 * (t311 * t158 * t181 - t53 * t181
     # * t162)))
      t308 = t149 + t24
      t309 = t275 * t1
      t215 = t215 * t26
      t317 = t35 * t1
      t325 = t18 * t181
      t363 = t259 * t160
      t364 = t74 * t8
      t365 = t61 * t160 * t250
      t366 = t8 * t162 * t181
      t367 = t366 - t330
      t368 = t20 * t367
      t369 = t207 * t66
      t183 = t365 * (t61 * (t25 * (t199 + t146) - t318 * t74) + t368 + t
     #369 * (t81 * t8 - t1) + t366 * (t49 - t35 - t37) + t329 * t2 * (-t
     #183 + t182)) + t365 * (t25 * (t61 * (t130 * t26 + t262) + t364 * (
     #-t51 - t19 - t20 - t38 - t49 - t50) * t181) + t338 * (t8 * (t154 *
     # (-t37 + t53) - t357) - t199 + t265 * t81) + t353 * (t364 * t334 +
     # t199 - t351) + (-t310 * t9 + t19 + t50 + t51) * t8 * t181 * t162 
     #+ t330 * (-t51 - t19 - t49 - t50)) + t363 * (t24 * (t1 * t308 + t8
     # * (-t190 * t203 + t181 * (t36 + t35) * t25)) + t32 * t181 * ((t25
     # * t323 + t74 * (t36 - t18)) * t61 * t8 + t309) + t325 * (t328 + t
     #215 - t317))
      t172 = t176 * (-t271 * t158 + t362 * t345 + t172)
      t203 = t30 + t29
      t262 = t44 * t53 * (t25 * (t61 * (-t329 * (t49 + t38) - t29) + t74
     # * (t354 - t32 + t35) * t181) + t325 * t74 * t267 + t329 * t203) +
     # t362 * (t1 * (t24 * (t282 * (t51 + t19 + t20 - t37 + t348 + t49 +
     # t50) - t311 * (t51 + t19 + t20 + t36 + t50)) - t44 * t27 + t321 *
     # (t74 * (t304 * t61 + t24) - t29 - t30)) + t313 * (-t260 * t211 * 
     #t250 - t1 * t203 + t81 * (t8 * (t30 + t28) - t223)))
      t265 = t154 * t248
      t304 = (-t37 - t35 + t36) * t8
      t325 = t336 * t159
      t328 = t32 - t18
      t330 = -t310 + t1
      t334 = t154 * t36
      t345 = t108 * t61
      t190 = t325 * (t1 * (-t358 + t50) + t8 * (t56 * (t149 * t328 + t36
     # * t24) + t29 - t334 + t190 * t181) + t49 * t330 + t51 * t330 + t2
     #0 * t330 + t19 * t330 + t345 * t315) + t325 * (t8 * (-t329 * t50 +
     # t154 * (t329 * t266 - t35)) + t349 * (t51 + t19 + t20 + t35 + t37
     # + t265 + t49 + t50) + t353 * (t251 + t177 + t304))
      t348 = t149 + t24
      t354 = t44 * t159 * t181 * t250 * t361
      t140 = -t126 * t356 * t158 * (t61 * (t35 + t53) + t24) - t140 * (t
     #307 + t161 * (t25 * (t56 * (t24 * (t214 * t26 - t199) + t53 * (-t1
     #0 * t44 * t181 + t211 * t66 * (t38 + t20)) + t249 * t61 * t267) + 
     #t49 * t61 * (t211 * (t282 - t311) + t44)) + t338 * (t74 * t290 * t
     #8 + t322) + t258 * t61 * (t310 * t267 + t346) + t320 * (t61 * (t37
     # * t74 + t319) - t29)) + t260 * (t272 * (t1 * (t154 * t298 + t29) 
     #+ t81 * (t25 * (t266 * t61 * t1 - t211) + t18 * (-t335 - t1))) + t
     #154 * (t269 * (t1 * (-t10 * t81 - t35 + t36) - t146) + t281 + t211
     # * (t19 * t53 * t259 + (t35 - t19 - t20) * t25 * t250)))) + t141 *
     # t183 + t283 * t173 * t74 * t348 + t284 * t360 * t160 * t359 + t28
     #5 * t354 - t286 * t190 + t287 * t336 * t53 * t159 * t332 + t57 * t
     #262 - t172
      t141 = t44 * t192
      t172 = t1 * t24
      t183 = t172 * t136 * t98
      t171 = t176 * t171
      t190 = (t1 * (t24 * (t4 * (-t51 - t19 - t20 + t37 - t49 - t50) + t
     #305) - t306) + t223 * (-t37 + t18 - t33 + t34 - t35) + t18 * (-t37
     # * t1 + t230)) * t98
      t214 = t34 * (t51 + t50 + t49 - t32) + t32 * (-t36 + t18)
      t230 = t51 + t50 + t49 + t38 + t36 + t20 + t19
      t249 = t36 * t18
      t262 = t10 * t125
      t267 = t98 * t41
      t269 = -t33 + t18 + t32
      t272 = t38 + t36
      t281 = t212 - t211
      t237 = t237 * t281
      t282 = t49 * t162 * t281
      t236 = t236 * (t209 - t204)
      t284 = t148 * t5
      t285 = t29 * t44 * t269
      t184 = t267 * (t25 * (t25 * (t212 * t19 + t238) + t4 * (-t212 * t2
     #03 - t316 * t24) + t231 + t179 * t125 * t4 * (t38 + t19 + t20) + t
     #105 * t228) + t18 * (-t2 * t184 * t25 + t100 * (t262 - t29)) + t22
     #3 * (t25 * (t49 + t20 + t36) + t4 * (-t262 + t29 + t30))) + t131 *
     # (t32 * (t51 + t50 + t19) + t249 + t241) * t98 + t301 * t191 * t21
     #4 + t44 * (t24 * (t230 * t25 - t244) + t30 * (t32 + t18) + t34 * (
     #t30 + t29)) + t191 * t4 * (t8 * (t32 * t33 - t207) - t284 + t263 *
     # t281) + t267 * (t25 * (t24 * (t99 * (t37 - t19) - t194) + t44 * t
     #272 * t5) + t245 + t237 + t282 + t236 + t322 * t32) + t285
      t203 = -t133 - t132
      t231 = t26 * t256
      t244 = t44 * t6
      t262 = t180 - t148
      t286 = t207 * t100
      t287 = t261 * t100
      t290 = t49 * t4
      t298 = t51 * t4
      t305 = t50 * t4
      t306 = t106 * t6
      t307 = t12 * t100
      t316 = t252 * t118
      t248 = t34 * t248
      t218 = t218 * t4
      t189 = t252 * t189
      t121 = t191 * (t8 * (t25 * (t218 - t248) - t189) + t196 + t197 + t
     #18 * t104 * t210) + t121 * (t118 * (t4 * (t8 * (t26 * (t7 * (t34 -
     # t33) - t106 * t5) + t261 * t25) + t287 + t290 * t213 + t298 * t21
     #3 + t305 * t213) + t18 * (t8 * (t203 * t7 + t306) - t307)) + t316 
     #* (t4 * (-t51 - t19 - t20 - t35 - t49 - t50) - t113))
      t129 = t125 * (t1 * (t5 * (t162 * t266 - t295) + t270 * t145 * t98
     # * t240) + t18 * (t32 * t227 + t148 - t180) - t233 * t32) + t125 *
     # (t222 * (t49 * t213 + t50 * t213 + t51 * t213 + t217 + t212 * (t2
     #67 * t12 + t37) + t211 * (t38 + t19)) + t221 + t32 * (t25 * (t1 * 
     #t129 - t243) - t216))
      t196 = t32 - t18
      t197 = (0.9q1 / 0.4q1)
      t210 = -(0.5q1 / 0.8q1)
      t216 = (0.5q1 / 0.16q2)
      t217 = 6 * t25
      t221 = 4 * t18
      t201 = t202 + t201
      t202 = t201 * t162
      t207 = t207 * t8
      t102 = (t25 * (t25 * (t175 + t238) + t258 * t281 + t105 * t281 + t
     #235 * t281 - t193 * t201 * t4) + t236 + t237 + t282 + t32 * t193 *
     # (t209 + t204)) * t98 * t40 + t202 * t193 + t299 * t145 * (t18 * t
     #281 + t252 * t24) + t301 * t299 * t214 + t110 * (t241 + t36 * (-t3
     #3 + t18) + t32 * (t51 + t19 + t49 + t50)) * t98 + t44 * (t25 * (t2
     #4 * (t51 + t19 + t20 + t36 + t49 + t50) + t29 * t5) + t30 * t269) 
     #+ t299 * t4 * (t8 * (t34 * (t38 + t19 + t20) + t302) - t284 - t207
     # - t10 * t223) + t285 + (t25 * (t25 * (t297 + t296) + t107 * (-t21
     #2 - t211)) + t245 + t18 * t2 * (t102 * t26 - t244) + t223 * (t26 *
     # (t186 + t106) + t105)) * t98 * t40
      t107 = t32 * (t51 + t19 + t49 + t50)
      t106 = t186 - t106
      t175 = t4 * t40 * t98
      t186 = t145 * t240
      t10 = t152 * (t32 * (-t44 * t114 + t232 - t301 * (t29 + t28)) + t4
     #4 * (t186 * t40 * t98 + t10 * (t33 - t18))) + t152 * (t1 * (t34 * 
     #(t12 * (t175 - t25) + t220) - t189) + t18 * (t1 * (t147 + t291) - 
     #t180) + t294 * t213 + t292 * t213 + t293 * t213 + t274 * t4 + t219
     # * t222)
      t114 = t299 * (t8 * (t25 * (-t218 + t248) + t276 * t226) + t195 + 
     #t224 * t20 + t224 * t19) + t151 * (t118 * (t4 * (t8 * (t25 * (t7 *
     # (-t247 - t246) - t264) + t231 * t4) - t287 + t290 * t228 + t298 *
     # t228 + t305 * t228) + t18 * (t8 * (t7 * (t133 + t132) - t306) + t
     #307)) + t316 * (t4 * (t51 + t19 + t20 + t35 + t49 + t50) + t113))
      t118 = (t138 - t130) * t25 * t26
      t132 = t44 * t230
      t133 = t56 * (-t18 * t61 - t24) + t154
      t151 = -t49 - t38 - t20 - t19
      t180 = t8 * t160 * t181 * t250
      t189 = -t176 - t180
      t195 = (t1 * t2 + t329 * (t138 + t130)) * t26
      t103 = t154 * (-t1 * t103 - t310 * t144)
      t144 = t289 * t25
      t201 = t29 * t1
      t212 = -t51 - t50 - t36
      t213 = t26 * (t255 + t254)
      t214 = t29 * t33
      t218 = t29 * t25
      t207 = t207 * t336
      t219 = t250 * t61
      t220 = -t271 * t211 + t318 * t66
      t224 = t141 * t4
      t226 = t25 * (-t357 * t1 * t136 + t158 * t220 + t224)
      t227 = t154 * t266 - t19 - t20 + t37 - t49 - t50 - t51
      t230 = t242 * t61 - t24
      t232 = (-t30 - t28) * t56 * t8 + t44
      t233 = t34 * t24
      t180 = t25 * (t1 * (t25 * (t37 * (-t149 - t24) + t233) + t357 * t2
     #27 - t185 * t308) + t18 * (t149 * t232 + t176 + t180) + t223 * (t7
     #4 * t230 + t25 * t308))
      t215 = t229 * t342 + t18 * (t154 * (t304 + t177) + t215 - t317) + 
     #t32 * (t327 - t351 + t309)
      t138 = t1 * (t329 * t136 + t334) + t368 + t199 * t61 * (t313 + t25
     #) + t138 * (-t162 * t6 * t181 + t154 * t26 - t329 * t26) + t50 * t
     #367 + t49 * t367 + t51 * t367 + t19 * t367 - t369 * t1 - t338 * t1
     #99 + t345 * (-t154 * t27 + t329 * (-t314 + t27))
      t199 = t158 * t56
      t229 = t260 * t250 * t215 + t199 * t138
      t236 = t353 * (t251 + t324 + t177) + t340 + t349 * (t51 + t19 + t2
     #0 + t35 + t265 + t49 + t50) + t329 * (t8 * (-t51 - t19 - t20 + t34
     #1 - t49 - t50) - t177)
      t104 = t61 * (t27 * (t328 * t56 - t25) + t209 * (-t329 + t154)) * 
     #t104 + t146 + t352
      t146 = t310 * (t51 + t20 + t38 + t50)
      t154 = t177 * t61
      t177 = t235 * t330
      t237 = t288 * t330
      t234 = t321 * (t8 * (t25 * (t333 + t24) + t29 + t30) + t234 * t61 
     #+ t18 * (t335 + t1))
      t238 = t276 * t1 * t56 * t315
      t134 = t366 * t134
      t240 = t320 * t24
      t199 = t199 * (t355 * (t51 - t53 + t20 + t32 - t37 + t50) + t134 -
     # t240) + t161 * (t25 * (t154 * t268 - t146) + t177 + t237 + t234 +
     # t238)
      t241 = t176 * t53 * (t56 * (-(0.9q1 / 0.2q1) * t53 + t221) + t217)
      t81 = t126 * t229 + t197 * t355 * t158 * t359 + t210 * t325 * t8 *
     # t196 + t216 * t219 * t159 * t332 + t283 * t199 - t57 * (t365 * t1
     #04 + t365 * t236) + t1 * (t25 * (t25 * (t149 * t212 - t213) + t214
     # - t288 * t348 - t105 * t348 - t235 * t348 + t172 * (t74 * t9 + t2
     #04)) + t18 * (t74 * (t30 - t161) - t218)) + t25 * (t53 * (t250 * t
     #158 * t8 * t151 * t181 + t103 + t195) + t50 * t189 + t51 * t189) +
     # t1 * (t1 * (t156 * t27 - t202) + t24 * t160 * t250) + t32 * t81 *
     # (t260 * t272 * t250 * t8 - t132) + t18 * t81 * (t118 + t225 - t38
     # * (t260 * t8 * t250 + t44)) - t211 * t81 * t162 * t143 + t18 * (t
     #363 * t211 * t61 + t201 * t74 - t144) + t32 * t74 * (t161 * t8 * t
     #133 - t201) - t245 * t81 + t207 * t160 + t223 * (-t30 * t165 - t21
     #8) + 3 * t226 - 2 * t180 + t241
      t156 = t344 * t8
      t158 = -t176 - t156
      t159 = t25 * (-t312 * t1 * t136 + t167 * t220 + t224)
      t160 = t150 + t24
      t156 = t25 * (t1 * (t25 * (t37 * (-t150 - t24) + t233) + t312 * t2
     #27 - t185 * t160) + t18 * (t150 * t232 + t156 + t176) + t223 * (t1
     #60 * t25 + t58 * t230))
      t160 = t167 * t56
      t138 = t280 * t250 * t215 + t160 * t138
      t134 = t160 * (t355 * (t51 - t55 + t20 + t32 - t37 + t50) + t134 -
     # t240) + t170 * (t25 * (t154 * t178 - t146) + t177 + t234 + t237 +
     # t238)
      t146 = t176 * t55 * (t56 * (-(0.9q1 / 0.2q1) * t55 + t221) + t217)
      t154 = 0.1q1 / t15
      t16 = 0.1q1 / t16
      t160 = -0.1q1 / t86
      t161 = -0.1q1 / t73
      t42 = 0.1q1 / t42
      t176 = 0.1q1 / t153
      t177 = 0.1q1 / t27
      t21 = 0.1q1 / t21
      t14 = 0.1q1 / t14
      t178 = t117 ** 2
      t180 = t137 ** 2
      t189 = t176 ** 2
      t199 = t63 ** 2
      t215 = t23 * t137 * t16
      t220 = t9 * t180 * t154
      t224 = t23 * t25
      t95 = 0.1q1 / t95
      t96 = 0.1q1 / t96
      t226 = 0.1q1 / t40
      t227 = 0.1q1 / t65
      t86 = 0.1q1 / t86
      t73 = 0.1q1 / t73
      t70 = 0.1q1 / t70
      t94 = 0.1q1 / t94
      t67 = 0.1q1 / t67
      t229 = 0.1q1 / t41
      t230 = 0.1q1 / t85
      t69 = 0.1q1 / t69
      t232 = t17 * t46
      t233 = t3 * t47
      t234 = t233 + t232
      t237 = t60 ** 2
      t238 = t2 * t7
      t240 = t28 * t16
      t241 = t22 * t154
      t242 = t13 * t229
      t22 = t22 * t46
      t243 = t77 * t78
      t246 = t243 * t75
      t247 = t246 * t86
      t248 = t22 * t3
      t251 = t93 * t95
      t252 = t21 * t3
      t254 = t252 * t16 * t23
      t255 = t240 * t66
      t258 = t14 * t42
      t163 = t258 * (t255 * t67 * t63 * t154 * t71 * t23 * t3 * (t47 * t
     #237 * t62 * t227 * t73 * t70 * t68 - t247 * t230 * t69 * (t47 * t7
     #5 + t22) + t22 * t60 * t62 * t227 * t73 * t70 * t68) + t252 * (t24
     #2 * (-t241 * t2 + t240) * t226 * t46 * t3 + t177 * t25 * (t1 * (t2
     #20 * t63 * t176 * t177 - t215 * t63 * t176) + t238 * t154)) + t254
     # * (t251 * (t234 * t90 - t248) + (t112 * (-t233 - t232) + t248) * 
     #t94 * t39) * t98 * t96 * t31 + t64 * t23 * t117 * (t163 + t72 * (-
     #t34 - t32 + t33 - t18) * t25 + (-t34 + t32 - t18) * t4 * t239 * t2
     #7 ** 2) * t180 * t181 * t154 * t16 * t250 * t177 * t176 * t160 ** 
     #2 * t161 ** 2)
      t227 = t60 * t62
      t230 = t251 * t91
      t233 = t124 * t39 * t94
      t240 = t28 * t66
      t248 = t47 * t154
      t260 = t14 * t16
      t15 = t42 * t3 * (t260 * t23 * (t240 * t67 * t71 * (-t227 * t73 * 
     #t70 * t68 + t247 * t69) + t252 * (t233 - t230) * t98 * t96 * t31 *
     # t15) + t242 * t21 * (t25 * t46 * t16 + (t248 * t25 + t1) * t14 * 
     #t3) * t226 * t2)
      t247 = 0.1q1 / t1
      t85 = -0.1q1 / t85
      t65 = -0.1q1 / t65
      t261 = t48 * t16
      t264 = t261 + t248
      t265 = t128 * t47 + t22
      t153 = 0.1q1 / t153 ** 2
      t266 = t258 * t3
      t22 = t266 * (t66 * t63 * t16 * t154 * t71 * t23 * (t1 * (t22 * t1
     #17 * t137 * t160 * t161 * t63 * t189 - t47 * t178 * t137 * t160 * 
     #t161 * t63 * t189) - t72 * t128 * t139 * t85 * t65 * t153 * t265) 
     #+ t252 * (t247 * (t177 * t2 * (t25 * t264 + t46 * (-t16 * t26 - t2
     #41)) + (t154 * t265 + t16 * (t128 * t48 - t46 * (t182 * t63 - t26)
     #)) * t153 * t239 * t139 * t27) + t177 * (t117 * t264 + t46 * (-(t7
     #2 * t6 + t26) * t16 - t241)) * t189 * t199 * t137 * t1))
      t47 = t78 ** 2
      t63 = t78 * t47
      t65 = t226 ** 2
      t72 = t96 ** 2
      t85 = t96 * t72
      t128 = t68 ** 2
      t137 = t68 * t128
      t139 = t2 ** 2
      t153 = t77 ** 2
      t239 = t62 ** 2
      t247 = t13 ** 2
      t264 = t229 ** 2
      t265 = t67 ** 2
      t267 = t97 * t94
      t268 = t267 * t142 * t39
      t269 = t251 * t92
      t270 = t269 * t111
      t271 = t268 * t264
      t274 = t23 * t16
      t276 = t274 * t98
      t281 = t60 * t239
      t282 = t75 * t153
      t255 = t258 * (t255 * t23 * t154 * t67 * t265 * t71 * (-t281 * t15
     #5 * t73 * t70 * t128 + t282 * t164 * t86 * t47 * t69) + t252 * (t3
     #1 * (t7 * t154 * (t270 * t65 + t271) * t72 + t276 * (-t270 * t226 
     #+ t268 * t229) * t85) + t139 * t13 * t247 * t45 * t65 * t154 * t26
     #4 * (t226 + t229)))
      t268 = t97 ** 2
      t270 = t39 ** 2
      t284 = t92 ** 2
      t285 = t70 ** 2
      t287 = t95 ** 2
      t289 = t93 ** 2
      t290 = t31 ** 2
      t291 = t94 ** 2
      t292 = t60 * t155 * t61 * t128
      t293 = t47 * t69
      t294 = t293 * t86
      t295 = t229 * t94 * t39
      t296 = t90 * t91
      t297 = t268 * t112
      t298 = t252 * t154
      t34 = t258 * (t298 * t13 * (t139 * t13 * (t44 * t52 - 2 * t44 * (t
     #37 + t33 + t35)) * t264 * t65 + (t297 * t124 * t270 * t291 * t264 
     #- t296 * t284 * t287 * t65 * t289) * t96 * t290) + t274 * (t24 * (
     #t186 * t25 * t162 * t139 * t154 * t47 * t128 * t177 + t311 * t265 
     #* t71 * (t239 * (t73 * (-t292 * t285 + (-t60 * (t25 * (t34 - t33 -
     # t55) + t32 * t157 + t18 * (-t58 + t25)) * t61 + t155) * t128 * t7
     #0) - t292 * t70 * t73 ** 2) + t294 * t153 * ((-t164 * (t86 + t69) 
     #- t25 * (t34 - t33 - t53) - t32 * t165 - t18 * (-t74 + t25)) * t61
     # * t75 + t164)) * t154) + t252 * (-t2 * t45 * t247 * t264 * t65 + 
     #(t295 * (-t127 * t97 + t142 * (t267 - t7)) + t251 * (-t101 * t92 -
     # t111 * t7) * t226 + t92 * t111 * t287 * t226 * t93) * t98 * t72 *
     # t31)))
      t37 = 0.1q1 / t84
      t52 = 0.1q1 / t82
      t53 = 0.1q1 / t88
      t74 = 0.1q1 / t89
      t82 = t267 * t112 * t39
      t84 = t82 * t229
      t88 = t269 * t90
      t89 = t88 * t226
      t165 = t89 - t84
      t292 = t3 ** 2
      t246 = t246 * t37 * t52 * t69
      t227 = t227 * t70
      t302 = t227 * t68 * t53 * t74
      t304 = t96 * t31
      t305 = t44 * t226
      t306 = t242 * t25
      t307 = t242 * t11
      t48 = t42 * (-t307 * t162 * t46 * t139 * t154 * t16 * t78 * t68 * 
     #t226 + t304 * t260 * t48 * t21 * t165 * t292 + t154 * t3 * t14 * (
     #t306 * (t21 * (t130 + t232) + t54 * t46 * t16 * t68 * t78) * t226 
     #* t2 + (t44 * t65 * t229 * t21 + t305 * t264 * t21) * t247 * t139 
     #+ t338 * t98 * t16 * t46 * (t28 * (-t302 + t246) * t67 + t304 * (t
     #88 * t37 * t74 - t82 * t52 * t53))))
      t54 = t75 ** 2
      t82 = t69 ** 2
      t88 = t90 ** 2
      t130 = t112 ** 2
      t232 = t239 * t70
      t260 = t54 * t153
      t308 = t237 * t239
      t309 = t269 * t88 * t37
      t310 = t295 * t268
      t311 = t97 * t229
      t17 = t42 * (t154 * t16 * t25 * (t98 * t46 * (t237 * (t232 * t128 
     #* t53 * t74 - t3 * t62 * (t150 * t57 + t24) * t70 * t128 * t53 * t
     #74 * t14) + t293 * t54 * t77 * t37 * t52 * (t3 * (t149 * t57 + t24
     #) * t14 - t77)) * t67 * t24 + t224 * t71 * t14 * (-t308 * t73 * t2
     #85 * t137 + t260 * t86 * t63 * t82) * t67 * t192 - t307 * t305 * t
     #23 * t139 * t78 * t68 * t14) + t254 * t14 * (t311 * t112 * t124 * 
     #t270 * t291 - t296 * t92 * t287 * t226 * t289) * t96 * t290 + t304
     # * t154 * t46 * (t66 * (t310 * t130 * t52 * t53 + t309 * t226 * t7
     #4 * (t3 * t83 * t14 - t92) - t267 * t3 * t130 * t80 * t39 * t52 * 
     #t53 * t229 * t14) * t16 + t252 * t17 * t14 * (-t89 + t84)))
      t86 = t66 * t31
      t89 = t269 * t226
      t149 = t66 * t7
      t254 = t268 * t39 * t94 * t264
      t45 = t258 * (t298 * (-t238 * t45 * t247 * t264 * t65 + (t284 * (t
     #111 * t287 * t65 * t93 - t251 * t101 * t65) + t254 * (t142 * t94 -
     # t127)) * t72 * t31) + t274 * (t28 * t154 * t71 * (t60 * (-t149 * 
     #t62 * t155 * t73 * t70 * t128 + t232 * t155 * t73 * t137) + t294 *
     # t75 * t77 * t164 * (-t149 + t243)) * t265 + t304 * (t86 * t130 * 
     #t124 * t268 * t154 * t270 * t229 * t52 * t53 * t291 - t271 * t252 
     #* t96 - t89 * (t86 * t230 * t88 * t92 * t37 * t154 * t74 + t252 * 
     #t111 * t96 * t226)) + t154 * t98 * (t308 * t59 * t285 * t137 * t53
     # * t74 - t260 * t76 * t37 * t52 * t63 * t82) * t67 * t162 * t192))
      t73 = t74 ** 2
      t82 = t74 * t73
      t86 = t53 ** 2
      t88 = t52 ** 2
      t101 = t37 ** 2
      t127 = t90 * t79
      t149 = t127 * t101 * t226
      t155 = t149 * t61 * t73 * t93 * t287
      t10 = -t155 + t251 * (t90 * (t126 * (t152 * (t100 * (t33 * t136 + 
     #t188) - t286 - t206 + t200 + t301 * (t26 * (-t33 * t2 + t25 * t256
     #) + t35 * t28) + t18 * (t8 * (t106 * t26 - t185) - t109 + t194 + t
     #182 * t106) + t32 * (t99 * (-t35 + t20) - t194 + t100 * (-t36 + t2
     #0))) + t152 * (t8 * (t162 * (-t12 * t24 + t323 * t5) + t107 * t25)
     # + t100 * (t33 * t300 + t107 + t249) + t109 * t120)) - t197 * t152
     # * t44 * t4 * t273 + t210 * t208 * t116 * t196 - t216 * t116 * t4 
     #* (t5 * (t8 * (t347 * t4 - t113) - t122 - t123) - t277 + t278) - t
     #283 * t10 - t57 * t114 + 3 * t222 * (t152 * t18 * t228 - t183 * t4
     #0 + t141) - 2 * t25 * (t299 * t263 * t179 + t190 * t40 + t171) - t
     #102 + t110 * t33 * t98 * (t5 * (-(0.9q1 / 0.2q1) * t175 + t217) + 
     #t221)) * t61 + t79) * t73 * t226 * t101
      t5 = t10 * t284 + t310 * t86 * t88 * ((-t119 * t94 + t126 * (t125 
     #* (t8 * (t257 * t162 * t6 - t29 * t147) + t200 - t286 + t50 * t262
     # + t49 * t262 + t51 * t262 + t20 * t262 + t19 * t262 + t35 * t148 
     #+ t108 * (t198 * t27 - t28 * t9) + t32 * (-t109 - t194) + t18 * (-
     #t109 + t194)) + t125 * (t4 * (t1 * (t145 * t192 + t188) + t99 * (t
     #231 + t303)) + t18 * (t2 * (-t100 * t6 + t8 * t203) + t7 * (t135 *
     # t26 + t244)) + t32 * (t100 * t275 + t326))) - t283 * t129 + t57 *
     # t121 + 3 * t222 * (t174 * t228 - t183 * t41 + t141) + t216 * t253
     # * t4 * t279 - t197 * t125 * t44 * t4 * t273 + t210 * t208 * t253 
     #* t196 - 2 * t25 * (t263 * t191 * t179 + t190 * t41 + t171) - t184
     # + t131 * t33 * t98 * (t5 * (-(0.9q1 / 0.2q1) * t187 + t217) + t22
     #1)) * t61 * t112 + t119)
      t1 = t60 * (t126 * t138 + t197 * t280 * t44 * t359 + t210 * t350 *
     # t8 * t196 + t216 * t219 * t168 * t332 + t283 * t134 - t57 * (t339
     # * t104 + t339 * t236) + t1 * (t25 * (t25 * (t150 * t212 - t213) +
     # t214 - t288 * t337 - t105 * t337 - t235 * t337 + t172 * (t58 * t9
     # + t204)) + t18 * (-t218 + t58 * (-t170 + t30))) + t25 * (t55 * (t
     #250 * t167 * t8 * t151 * t181 + t103 + t195) + t50 * t158 + t51 * 
     #t158) + t1 * (t1 * (t166 * t27 - t202) + t331 * t24) + t343 * (t28
     #0 * t272 * t250 * t8 - t132) + t18 * t87 * (-t38 * (t280 * t8 * t2
     #50 + t44) + t118 + t225) - t211 * t87 * t162 * t143 + t18 * (t211 
     #* t169 * t61 * t259 + t201 * t58 - t144) + t32 * t58 * (t170 * t8 
     #* t133 - t201) - t245 * t87 + t207 * t169 + t223 * (-t30 * t157 - 
     #t218) + 3 * t159 - 2 * t156 + t146)
      t6 = t140 * t69
      t8 = t293 * t88 * t101
      t9 = t28 * t265
      t10 = t274 * t258 * t154
      t1 = t10 * (t186 * t25 * t139 * t173 * t43 * t247 * t47 * t128 * t
     #264 * t65 + t9 * (t239 * (-t60 * t205 * t61 * t128 * t86 * t73 * t
     #285 + (t1 * t61 + t205) * t73 * t86 * t128 * t70) + t8 * t153 * ((
     #-t6 + t81) * t61 * t75 + t140)) * t56 * t115 + t219 * t72 * t31 * 
     #t5 * t98)
      t5 = t251 * t284
      t12 = t7 * t98
      t18 = t281 * t205 * t70
      t19 = t72 * t31
      t6 = t10 * (t9 * t115 * (t6 * t282 * t101 * t88 * t63 + t18 * t137
     # * t86 * t73) - t240 * t7 * t265 * t115 * (t227 * t205 * t128 * t8
     #6 * t73 + t8 * t75 * t77 * t140) + t19 * (t90 * (-t12 * t89 * t79 
     #* t101 * t73 + t5 * t79 * t101 * t65 * t73) + t84 * t119 * t88 * t
     #86 * (-t12 + t311)) * t250 * t181)
      t7 = t267 * t39
      t10 = t304 * t252 * t42 * (t46 * t16 * t165 + (t89 * (t248 * t90 +
     # t91) - t7 * t229 * (t248 * t112 + t124)) * t14 * t3)
      t11 = t258 * (-t274 * t11 * t162 * t139 * t154 * t78 * t68 * t177 
     #+ (-t261 * t306 * t2 * t226 + t304 * (-t295 * (t112 * t80 * t16 + 
     #t241 * t97) + t251 * (t90 * t83 * t16 + t241 * t92) * t226) * t46)
     # * t21 * t292 + t242 * t252 * t226 * t2 * (t135 * t154 * t2 + t274
     # * t44))
      t12 = t276 * t258 * t66 * t154
      t7 = t12 * (t28 * (-t237 * t62 * t70 * t68 * t53 * t74 * t234 + t2
     #43 * t234 * t69 * t52 * t37 * t54) * t67 + t304 * (-t7 * t130 * t5
     #2 * t53 * t234 + t309 * t234 * t74))
      t20 = t295 * t297 * t119 * t88 * t86
      t5 = t12 * (t31 * (t66 * (t5 * t149 * t73 - t20) * t85 + t66 * (t1
     #27 * (-t37 * t101 * t226 * t73 * t93 - t226 * t101 * t82 * t93) * 
     #t95 * t284 - t20 * (t52 + t53)) * t72) + t9 * t98 * (t18 * (t73 * 
     #t53 * t86 + t86 * (-t67 * t73 + t82)) * t128 + t8 * t282 * t140 * 
     #(t37 + t52 + t67)))
      ret = -12 * t266 * t64 * t215 * t66 * t117 * t176 * t160 * t161 * 
     #t71 - 16384 * t5 + 320 * t276 * t266 * t66 * (t28 * (-t246 * t76 +
     # t302 * t59) * t67 + t304 * (t233 * t97 * t112 * t52 * t53 - t230 
     #* t90 * t92 * t37 * t74)) + 1024 * t255 + 8192 * t6 + 64 * t48 + 9
     #6 * t10 + 256 * t34 + 4 * t22 + 48 * t15 - 16 * t163 - 8 * t14 * t
     #177 * t42 * (t66 * t193 * t23 * t178 * t180 * t199 * t154 * t16 * 
     #t189 * t160 * t161 + (t2 * (-t224 * t16 + (-t209 * t177 + t4) * t1
     #54 * t2) + (-t220 * t199 * t189 * t177 + t215 * t199 * t189) * t11
     #7 * t193) * t21 * t3) - 192 * t7 - 128 * t17 - 32 * t11 - 4096 * t
     #1 - 2048 * t19 * t298 * t258 * ((t13 * t226 * t65 * t93 - t96 * t6
     #5 * t93) * t95 * t111 * t284 + t254 * t142 * (t242 + t96)) - 512 *
     # t45

      hjetmass_qqbgg_bubble_pmpm_s123 = -ret/16q0*(0,1q0)
      return

      end function

      complex*32 function hjetmass_qqbgg_bubble_pmpm_s124
     &     (i1,i2,i3,i4,za,zb,mt,p,flip)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt
          real*16 p(mxpart,4)
          complex*32 alpha
          logical flip

      alpha = (za(i1,i2)*zb(i2,i1) + za(i1,i4)*zb(i4,i1) +
     & za(i2,i4)*zb(i4,i2))/(za(i1,i2)*zb(i2,i1) + za(i1,i4)*zb(i4,i1))

      p(5,:) = real(alpha)*p(i1,:)
      p(6,:) = (1q0-real(alpha))*p(i1,:) + p(i2,:) + p(i4,:)
      if (flip .eqv. .true.) then
          call spinoru_qp(6,p,zb,za)
      else
          call spinoru_qp(6,p,za,zb)
      end if

      t1 = za(6, i1)
      t2 = za(i2, i4)
      t3 = zb(i2, 5)
      t4 = zb(i3, i1)
      t5 = za(6, i3)
      t6 = za(i1, i2)
      t7 = za(i1, i3)
      t8 = zb(i2, i1)
      t9 = zb(i3, 5)
      t10 = za(i1, i4)
      t11 = zb(i4, i1)
      t12 = t6 * t8
      t13 = t10 * t11
      t14 = t2 * zb(i4, i2)
      t15 = t12 + t13 + t14
      t16 = zb(i4, 5)
      t17 = za(6, i2)
      t18 = za(5, i2)
      t19 = zb(i2, 6)
      t20 = za(5, i3)
      t21 = zb(i3, 6)
      t22 = za(5, i4)
      t23 = zb(i4, 6)
      t24 = t18 * t19
      t25 = t20 * t21
      t26 = t22 * t23
      t27 = t26 + t24 + t25
      t28 = zb(i1, 6)
      t29 = za(6, i4)
      t30 = t17 * t3
      t31 = t5 * t9
      t32 = t29 * t16
      t33 = t31 + t32 + t30
      t34 = t1 * t28
      t35 = t18 * t3
      t36 = t17 * t19
      t37 = t20 * t9
      t38 = t5 * t21
      t39 = t22 * t16
      t40 = t29 * t23
      t41 = -t39 + t40 - t37 + t38 + t34 - t35 + t36
      t41 = t41 ** 2
      t42 = 4 * t33 * t27 + t41
      t42 = sqrt(t42)
      t43 = t42 + t39 - t40 + t37 - t38 - t34 + t35 - t36
      t42 = t42 - t39 + t40 - t37 + t38 + t34 - t35 + t36
      t44 = 0.1q1 / t27
      t45 = (0.1q1 / 0.2q1)
      t46 = t45 * t22
      t47 = -t46 * t43 * t44 - t29
      t48 = t45 * t21
      t49 = -t48 * t43 * t44 + t9
      t41 = 4 * t33 * t27 + t41
      t41 = sqrt(t41)
      t50 = t41 + t39 - t40 + t37 - t38 - t34 + t35 - t36
      t50 = 0.1q1 / t50
      t51 = 2
      t52 = t45 * t43 * t44
      t53 = t51 * t33
      t54 = -t53 * t50 - t52
      t55 = t44 * (t43 + t42)
      t46 = t46 * t42 * t44 - t29
      t48 = t48 * t42 * t44 + t9
      t41 = t41 - t39 + t40 - t37 + t38 + t34 - t35 + t36
      t41 = 0.1q1 / t41
      t56 = t45 * t42 * t44
      t53 = t53 * t41 + t56
      t57 = 0.1q1 / t20
      t58 = t18 * t5
      t59 = -t58 * t57 + t17
      t60 = 0.1q1 / t21
      t61 = t9 * t60
      t62 = t5 * t57
      t63 = t61 + t62
      t64 = 0.1q1 / t23
      t65 = t16 * t64
      t66 = t65 + t62
      t67 = t65 * t18 + t17
      t68 = 0.1q1 / t22
      t69 = t29 * t68
      t70 = t65 + t69
      t71 = -t65 + t52
      t72 = -t65 - t56
      t73 = zb(i4, i3)
      t74 = t61 - t65
      t75 = -t69 * t18 + t17
      t76 = t69 + t52
      t77 = t69 - t56
      t78 = zb(i3, i2)
      t79 = -t52 * t18 - t17
      t80 = t56 * t18 - t17
      t81 = t5 * t59
      t82 = t20 * t17
      t58 = -t58 + t82
      t83 = t64 ** 2
      t84 = t16 ** 2
      t85 = t22 * t5
      t86 = -t85 * t57 + t29
      t87 = t62 - t69
      t88 = -t65 * t21 + t9
      t89 = -t52 * t19 + t3
      t90 = t56 * t19 + t3
      t91 = 0.1q1 / t27
      t92 = t7 * t4
      t93 = za(i3, i4) * t73
      t94 = za(i2, i3) * t78
      t95 = t26 + t24 + t25
      t96 = t91 ** 2
      t97 = t91 * t96
      t44 = (t40 + t92 + t93 + t94 + t12 + t13 + t14 + t34 + t36 + t38) 
     #* t44
      t98 = t95 * t96
      t99 = t98 * t43 + t44
      t100 = t39 - t92 - t93 - t94 - t12 - t13 - t14 + t35 + t37
      t101 = t43 ** 2
      t102 = t96 * t101
      t103 = t97 * t43 * t101
      t104 = (t39 + t35 + t37) * t91
      t105 = t104 * t43 + t30 + t31 + t32
      t106 = (0.1q1 / 0.4q1)
      t107 = (0.1q1 / 0.8q1)
      t40 = t40 + t92 + t93 + t94 + t12 + t13 + t14 + t34 + t36 + t38
      t108 = t40 * t17
      t109 = t18 ** 2 * t3
      t110 = (-t12 - t13 - t14 - t92 - t93 - t94 + t39 + t37) * t91
      t111 = t42 ** 2
      t112 = t98 * t111
      t113 = t95 * t91
      t114 = t113 * t42 + t35 + t37 + t39
      t115 = (0.3q1 / 0.4q1)
      t116 = t109 * t42 * t91 + t110 * t42 * t18 + t112 * t115 * t18 - t
     #51 * t17 * t114 + t108
      t117 = -t12 - t13 - t14 - t92 - t93 - t94 + t37 + t35
      t40 = t40 * t29
      t118 = t117 * t91
      t119 = t22 ** 2 * t16
      t26 = -t26 - t24 - t25
      t44 = t26 * t96 * t42 + t44
      t120 = t96 * t111
      t97 = t97 * t42 * t111
      t32 = t104 * t42 - t30 - t31 - t32
      t104 = t106 * t120 * t18 * t100 + t107 * t97 * t18 * t95 + t45 * t
     #42 * t17 * t44 - t17 * t32
      t98 = t98 * t101
      t26 = -t109 * t43 * t91 - t110 * t43 * t18 + t115 * t98 * t18 - t5
     #1 * t17 * (t26 * t91 * t43 + t35 + t37 + t39) + t108
      t101 = -t106 * t102 * t18 * (-t39 + t92 + t93 + t94 + t12 + t13 + 
     #t14 - t35 - t37) - t107 * t103 * t18 * t95 - t45 * t43 * t17 * t99
     # + t17 * t105
      t32 = t106 * t120 * t22 * t100 + t107 * t97 * t22 * t95 + t45 * t4
     #2 * t29 * t44 - t29 * t32
      t15 = 0.1q1 / t15
      t6 = 0.1q1 / t6
      t44 = 0.1q1 / t55
      t55 = 0.1q1 / t42
      t97 = 0.1q1 / t43
      t108 = 0.1q1 / t8
      t109 = 0.1q1 / t11
      t54 = 0.1q1 / t54
      t53 = 0.1q1 / t53
      t7 = 0.1q1 / t7
      t110 = t47 * t108
      t111 = t79 * t109
      t120 = t10 * t6
      t121 = t46 * t108
      t122 = t80 * t109
      t123 = t53 * t41
      t124 = t123 * t55
      t125 = t124 * t48
      t126 = t50 * t54
      t127 = t7 * t44
      t128 = -0.1q1 / t71
      t129 = -0.1q1 / t72
      t130 = t120 * t108
      t131 = t130 + t109
      t132 = t2 ** 2
      t133 = t78 * t108
      t134 = t123 * t48
      t135 = t29 * t108
      t136 = t17 * t109
      t137 = t97 * t55
      t138 = 0.1q1 / t63
      t72 = 0.1q1 / t72
      t71 = 0.1q1 / t71
      t139 = 0.1q1 / t5
      t140 = 0.1q1 / t66
      t141 = 0.1q1 / t70
      t66 = 0.1q1 / t66
      t142 = t46 * t73
      t143 = t78 * t80
      t144 = t126 * (-t47 * t73 + t78 * t79)
      t145 = t123 * (t143 - t142) + t144
      t146 = t135 + t136
      t147 = t108 * t5 * t86 + t109 * t81
      t148 = t138 ** 2
      t149 = t57 ** 2
      t150 = t148 * t140
      t69 = t69 * t141
      t151 = t91 * t44 * t33
      t152 = t34 * t64 * t6
      t153 = t137 * t27
      t154 = t109 * t108
      t155 = t154 * t28
      t12 = t7 * (t155 * (t152 * (t69 * (t88 * (-t51 * t65 * t17 * (t65 
     #* (t25 + t24) - t35 - t37) + t17 * (t31 + t30) + t117 * t18 * t83 
     #* t84 - t65 * (t12 + t13 + t14 + t92 + t93 + t94 + t38 + t34 + t36
     #) * t17 - t18 * (t25 + t24) * t64 * t83 * t16 * t84) * t72 ** 2 * 
     #t71 ** 2 * t96 - (t18 * t20 * t84 * t83 + t82 * t65 * t51 + t17 * 
     #t5) * t149 * t66 ** 2) + (-t140 ** 2 * t149 * t60 * t138 - t150 * 
     #t149 * t60) * t81 * t9) + t151 * t4 * (t123 * (t90 * t64 * t129 + 
     #t120) * t80 + t126 * t79 * (t89 * t64 * t128 + t120))) + (t151 * t
     #28 * (t109 * t145 + t130 * t145) + t34 * (-t6 * t147 * t60 * t139 
     #** 2 * t138 - t6 * t57 * t147 * t60 * t139 * t148) * t9 + t153 * t
     #4 * (t120 * t146 * t16 + t146 * t3)) * t15 * t2)
      t13 = 0.1q1 / t77
      t14 = 0.1q1 / t87
      t18 = 0.1q1 / t9
      t24 = 0.1q1 / t76
      t25 = t28 ** 2
      t30 = t60 ** 2
      t31 = t141 ** 2
      t36 = 0.1q1 / t70 ** 2
      t38 = 0.1q1 / t63 ** 2
      t63 = t68 ** 2
      t18 = t4 * t18
      t70 = t133 * t28
      t81 = t16 * t109
      t82 = t1 * t2
      t84 = t59 * t109
      t92 = t84 * t28 * t57
      t93 = t1 ** 2 * t2 * t4 * t6
      t94 = t57 * t66
      t117 = t82 * t6
      t3 = t64 * (t94 * t120 * t29 * t28 * t67 ** 2 * t63 * t31 * t109 -
     # t94 * t69 * t28 * t67 * t109) + t83 * (-t81 * t120 * t75 ** 2 * t
     #28 * t57 * t68 * t14 * t36 + t81 * t135 * t1 * t4 * t25 * t36 * t1
     #3 * t6 * t91 * t24 * (t10 * t75 - t82) * t63) + t9 * (t92 * t138 *
     # t60 * t140 * t64 + t82 * t30 * t6 * (t155 * t4 * t57 + (t70 * t57
     # + ((t62 * t23 + t16) * t109 - (t62 * t19 + t3) * t108) * t139 * t
     #4) * t15 * t2) * t148) + t69 * t155 * t64 * t91 * (t4 * t67 * (-t6
     #5 * t19 + t3) + t65 * (-t67 * t78 + t73 * (t65 * t22 + t29) + t120
     # * t1 * t4 * t67 * t141 * t68 - t93 * t141 * t68) * t28) * t71 * t
     #72 + t117 * (-t155 * t5 * t4 * t149 * t38 * t60 + ((-t70 * t60 + t
     #18 * ((-t61 * t23 + t16) * t109 - (-t61 * t19 + t3) * t108)) * t14
     #9 * t38 * t5 + t18 * t139 * (t108 * t3 - t81)) * t15 * t2)
      t13 = 0.1q1 / t74
      t18 = 0.1q1 / t74
      t19 = 0.1q1 / t87
      t24 = t15 * t2
      t3 = t7 * t3 + t7 * t109 * t28 * (t6 * (t57 * (t64 * (-t82 * t29 *
     # t67 * t31 * t66 * t63 + t61 * t59 * t138 * t140 * t19 * (t10 * t5
     #9 - t82) * t68) + t82 * t16 * (t135 * t28 * t73 * t13 * t141 * t66
     # * t60 + t75 * t36 * t14) * t68 * t83 + t132 * t1 * t9 * t73 * t14
     #8 * t15 * t30) + t82 * t60 * t73 * t5 * (-t24 * t38 + t61 * (-t18 
     #* t38 + t150) * t108 * t64 * t28) * t149) + t69 * t34 * t4 * t16 *
     # t108 * t83 * t91 * t71 * t72)
      t5 = -0.1q1 / t77
      t13 = -0.1q1 / t76
      t14 = t9 ** 2
      t18 = t123 * t42
      t19 = t126 * t43 * t128
      t36 = t18 * t129
      t38 = t88 * t72 * t71 * t91
      t59 = t6 * t1
      t62 = t155 * t64
      t66 = t62 * t96
      t70 = t7 * t28
      t5 = t70 * (t66 * (t36 * (-t143 + t142) + t93 * t68 * (-t18 * t5 *
     # t129 + t19 * t13) + (-t36 * (t120 * t80 * t5 * t68 + 1) + t19 * (
     #t120 * t79 * t13 * t68 + 1)) * t4 * t1 + t144 * t128 * t43) * t44 
     #* t33 + t59 * (t24 * ((t108 * t86 + t84) * t30 * t148 * t14 + t136
     # + t135) * t139 + t108 * t64 * (-t69 * t2 * t88 * t91 * t71 * t72 
     #- t28 * t67 * t109 * (t38 + t94) * t63 * t31 * t29 ** 2 + t92 * t1
     #50 * t14 * t30)))
      t13 = t34 * t91
      t14 = t13 * t108
      t18 = (-t110 - t111) * t97
      t30 = (t121 + t122) * t55
      t31 = t13 * t109
      t63 = t134 * t129
      t68 = t126 * t49
      t74 = t68 * t128
      t34 = t34 * t6
      t4 = t7 * (t34 * t9 * (t62 * t58 * t149 * t138 * t60 * t140 + t24 
     #* (-t153 * t146 + (t58 * t109 - (-t20 * t29 + t85) * t108) * t60 *
     # t57 * t139 * t138)) + t2 * (t152 * t91 * t108 * (t74 + t63) + (t1
     #26 * (t18 * t89 + t31) + t120 * (t123 * (t30 * (t56 * t23 + t16) +
     # t14) + t126 * (t18 * (-t52 * t23 + t16) + t14)) + t123 * (t30 * t
     #90 + t31)) * t15 * t4) * t44 * t33)
      t14 = t38 * t69 * t67
      t10 = t70 * t64 * (t14 * t109 + (t14 * t10 + t82 * (-t61 * t138 * 
     #t140 * t57 + t69 * (-t65 * t28 * t73 * t72 * t71 * t91 * t109 + t9
     #4))) * t108 * t6)
      t14 = t104 * t109 + t108 * t32
      t16 = (t100 * t102 * t106 * t22 - t103 * t107 * t22 * t95 - t45 * 
     #t43 * t29 * t99 + t29 * t105) * t108 + t101 * t109
      t20 = t44 ** 2
      t23 = t97 ** 2
      t31 = t55 ** 2
      t38 = t54 ** 2
      t22 = t34 * t24 * t7 * (t9 * t27 ** 2 * t23 * t31 * (t108 * t29 * 
     #t33 + t109 * t17 * t33) + (t124 * (t14 * t21 + t48 * (t108 * (t115
     # * t112 * t22 - t51 * t29 * t114 + t118 * t42 * t22 + t119 * t42 *
     # t91 - t32 * t53 + t40) + t109 * (-t104 * t53 + t116))) + (t54 * (
     #t16 * t21 + t49 * (t26 * t109 + (t115 * t98 * t22 - t118 * t43 * t
     #22 - t119 * t43 * t91 + t51 * t29 * (t113 * t43 - t35 - t37 - t39)
     # + t40) * t108)) - t49 * t16 * t38) * t97 * t50) * t91 * t20 * t33
     #)
      t32 = t41 ** 2
      t35 = t50 ** 2
      t37 = t53 ** 2
      t39 = t48 * t104
      t40 = t24 * t33
      t18 = t34 * t127 * t33 * (t40 * (t18 * t38 * t35 * t49 - t30 * t48
     # * t32 * t37) + t66 * (t41 * (t53 * (t129 * (t104 * t21 + t116 * t
     #48) - t39 * t129 ** 2) - t39 * t129 * t37) + t126 * t128 * (t101 *
     # (t49 * (t128 + t54) - t21) - t26 * t49)) * t44)
      t21 = t34 * t7 * t20 * t33 * (t66 * t44 * (t74 * t101 + t63 * t104
     #) + t24 * (t134 * t14 * t31 - t68 * t23 * t16))
      t23 = t59 * t154 * t127 * t25 * t33 ** 2 * t64 * t91 * (-t79 * t49
     # * t35 * t128 * t38 + t48 * t80 * t32 * t129 * t37)
      ret = 96 * t127 * t15 * t33 * t2 * (t126 * (-t111 * t8 - t47 + t12
     #0 * (-t110 * t11 - t79)) * t97 * t49 + t125 * (t122 * t8 + t46 + t
     #120 * (t121 * t11 + t80))) + 48 * t7 * ((t64 * (t126 * t131 * t128
     # * t79 * t49 + t134 * t80 * t129 * t131) + t6 * (-t133 * (t123 + t
     #126) + (-t123 - t126) * t109 * t73) * t15 * t132 * t1) * t91 * t44
     # * t33 * t28 + t137 * t15 * t27 * t9 * t2 * (t136 * t8 + t29 + t12
     #0 * (t135 * t11 + t17))) - 16 * t12 + 8 * t5 - 32 * t4 - 12 * t10 
     #+ 24 * t117 * t154 * t127 * t25 * t33 * t73 * t64 * t96 * (-t36 + 
     #t19) - 256 * t22 - 128 * t18 + 512 * t21 + 1024 * t40 * t13 * t7 *
     # t44 * t20 * t6 * (-t68 * t97 * t16 + t125 * t14) + 64 * t23 + 4 *
     # t3

      hjetmass_qqbgg_bubble_pmpm_s124 = -ret/16q0*(0,1q0)
      return

      end function

      complex*32 function hjetmass_qqbgg_bubble_pmpm_s12
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt

      t1 = za(i2, i4)
      t2 = zb(i3, i2)
      t3 = za(i2, i3)
      t4 = zb(i3, i1)
      t5 = zb(i4, i1)
      t6 = t3 * t4
      t7 = t1 * t5
      t8 = t7 + t6
      t9 = za(i1, i3)
      t10 = za(i1, i4)
      t11 = zb(i4, i2)
      t12 = t10 * t11
      t13 = t9 * t2 + t12
      t14 = t9 * t4
      t15 = t3 * t2
      t16 = t10 * t5
      t17 = t1 * t11
      t18 = -t17 + t14 - t15 + t16
      t18 = t18 ** 2
      t19 = 4 * t8 * t13 + t18
      t19 = sqrt(t19)
      t20 = -t17 + t14 - t15 + t16 + t19
      t21 = 0.1q1 / t13
      t22 = (0.1q1 / 0.2q1)
      t23 = t22 * t10
      t24 = -t23 * t20 * t21 - t1
      t25 = za(i3, i4)
      t18 = 4 * t8 * t13 + t18
      t18 = sqrt(t18)
      t26 = -t17 + t14 - t15 + t16 + t18
      t27 = 0.1q1 / t9
      t28 = t22 * t20 * t21
      t29 = t3 * t27
      t30 = -t28 - t29
      t26 = 0.1q1 / t26
      t31 = 2
      t32 = t31 * t8
      t33 = -t32 * t26 - t28
      t19 = -t17 + t14 - t15 + t16 - t19
      t34 = t21 * (-t20 + t19)
      t35 = zb(i4, i3)
      t23 = -t23 * t19 * t21 - t1
      t18 = -t17 + t14 - t15 + t16 - t18
      t36 = t22 * t19 * t21
      t37 = -t36 - t29
      t18 = 0.1q1 / t18
      t32 = -t32 * t18 - t36
      t38 = -t28 * t2 + t4
      t39 = t25 * t35
      t9 = t1 * t9
      t40 = t9 * t2
      t41 = t10 * (t15 + t39) - t40
      t42 = t3 * t10
      t9 = (-t9 + t42) * t4
      t43 = t39 * t1
      t44 = -t43 + t9
      t45 = -t22 * t20 * t21 * t41 + t44
      t46 = 0.1q1 / t11
      t47 = t5 * t46
      t28 = -t28 + t47
      t48 = -t36 * t2 + t4
      t36 = -t36 + t47
      t49 = t29 + t47
      t50 = 0.1q1 / t10
      t51 = t1 * t50 + t47
      t52 = t10 * t4
      t53 = t1 * t2
      t54 = t53 + t52
      t55 = 0.1q1 / t2
      t56 = t4 * t55 + t29
      t57 = -t42 * t27 + t1
      t58 = t15 * t27 + t4
      t59 = t16 * t46 + t1
      t60 = t2 * t5
      t61 = t60 * t46
      t62 = t4 - t61
      t13 = 0.1q1 / t13
      t63 = t14 - t15
      t64 = t10 ** 2
      t65 = t1 ** 2
      t66 = t65 * t2 * t11
      t5 = t64 * t4 * t5
      t67 = t53 * t63
      t63 = t52 * t63
      t68 = t60 * t1
      t69 = (t17 + t39) * t4
      t70 = t1 * (t68 + t69)
      t71 = (t5 - t66 + t67 + t63) * t13
      t72 = t13 ** 2
      t73 = t39 * t21 * (t53 - t52)
      t74 = t2 * (t10 * (t17 + t14) + t40) * t72
      t11 = t4 * t11
      t75 = t39 * t2
      t76 = t10 * (t11 + t60) + t75
      t77 = (0.1q1 / 0.4q1)
      t78 = t31 * t4 * (t10 * (t7 + t6) + t15 * t1)
      t79 = -t71 * t20 + t22 * t20 * (t74 * t20 + t73) + t77 * t72 * t20
     # ** 2 * t10 * t76 - t70 - t78
      t70 = -t71 * t19 + t22 * t19 * (t74 * t19 + t73) + t77 * t72 * t19
     # ** 2 * t10 * t76 - t70 - t78
      t14 = t2 * (t10 * (t17 + t14) + t40) * t13
      t11 = t11 + t60
      t17 = t2 * t19 * t13
      t15 = t10 * (t15 + t39) - t40
      t40 = t27 ** 2
      t3 = t3 ** 2
      t2 = t2 * t20 * t13
      t60 = 0.1q1 / t49
      t35 = 0.1q1 / t35
      t71 = -0.1q1 / t36
      t51 = 0.1q1 / t51
      t73 = -0.1q1 / t30
      t25 = 0.1q1 / t25
      t74 = -0.1q1 / t37
      t76 = -0.1q1 / t28
      t49 = 0.1q1 / t49
      t56 = 0.1q1 / t56
      t77 = t4 ** 2
      t78 = t62 ** 2
      t80 = t51 * t50
      t81 = t80 * t59
      t82 = t58 * t55
      t83 = t77 * t60
      t25 = t25 * t35
      t6 = t6 * t40 * t56
      t3 = t25 * (t7 * t80 * t46 ** 2 * (t54 * t27 * t49 + t78 * (t47 * 
     #t15 + t43 - t9) * t72 * t71 ** 2 * t76 ** 2) + t6 * t55 * (t54 * t
     #60 * t46 + t57 * (t31 * (-t16 * t1 * t4 + t29 * (t5 - t66) + t12 *
     # t53 * t3 * t40) + t1 * (-t68 - t69) + t29 * t39 * (t52 - t53) + t
     #10 * (t10 * t11 + t75) * t40 * t3) * t72 * t73 ** 2 * t74 ** 2))
      t7 = 0.1q1 / t34
      t9 = 0.1q1 / t36
      t12 = 0.1q1 / t30
      t16 = 0.1q1 / t28
      t28 = 0.1q1 / t33
      t30 = 0.1q1 / t37
      t32 = 0.1q1 / t32
      t33 = t48 ** 2
      t34 = t23 * t70
      t35 = t34 * t27 * t30
      t21 = (-t22 * t19 * t21 * t41 + t44) * t9
      t22 = t21 * t33 * t46 - t35
      t36 = t38 ** 2
      t37 = t7 ** 2
      t41 = t24 * t79 * t27 * t12
      t43 = t19 * t18
      t44 = t43 * t32
      t54 = t28 ** 2
      t68 = t26 ** 2
      t69 = t32 ** 2
      t75 = t24 ** 2
      t84 = t12 * t27
      t85 = t28 * t26
      t86 = t85 * t20
      t87 = t48 * t46
      t88 = t25 * t13
      t89 = t88 * t7 * t8
      t5 = t89 * (t8 * (t84 * t75 * t38 * t68 * t54 - t23 * t48 * t18 **
     # 2 * t69 * (t23 * t27 * t30 + t87 * t9) + t24 * t36 * t46 * t68 * 
     #t16 * t54) + (t43 * (t22 * t69 + t32 * (t27 * (t30 * (t10 * t70 + 
     #t23 * (-t31 * (t14 * t19 - t5 - t63 + t66 - t67) - t64 * t19 * t13
     # * t11 - t39 * (t10 * (t17 - t4) + t53))) - t34 * t30 ** 2) + t33 
     #* t9 * t46 * (t21 - t15))) + t86 * (t84 * (t24 * (-t31 * (t14 * t2
     #0 - t5 - t63 + t66 - t67) - t64 * t20 * t13 * t11 - t39 * (t10 * (
     #t2 - t4) + t53)) + t79 * (-t24 * (t12 + t28) + t10)) + (t45 * t16 
     #** 2 + t16 * (t28 * t45 - t15)) * t46 * t36)) * t72 * t7)
      t9 = t25 * t72
      t2 = t9 * t37 * t8 * (t85 * (t38 * t45 * t46 * t16 * (-t2 + t38) -
     # t41) + (t87 * t21 * (-t17 + t48) - t35) * t32 * t18)
      t10 = t57 * t4
      t11 = t53 + t52
      t9 = t9 * t27 * t7 * t8 * (-t86 * t24 * t12 * t11 + t44 * t11 * t3
     #0 * t23)
      ret = -8 * t25 * (t57 ** 2 * t77 * t58 * t27 * t56 ** 2 * t55 ** 2
     # * t13 * t73 * t74 + t83 * t27 * t46 * t55 * t56 * (t82 * t56 - 1)
     # * t57 + t46 * t27 * t4 * ((t82 + t29) * t56 * t60 + t51 * t49 * (
     #t50 * t59 + t47)) * t1 + t80 * t46 * (-t81 * t78 * t13 * t71 * t76
     # + (-t62 * (t81 + 1) + t61) * t49 * t27) * t65 + t83 * t42 * t40 *
     # t56 * t55 * t46) + 16 * t3 + 256 * t25 * t13 * t72 * t7 * t37 * t
     #8 * (t44 * t22 + (-t36 * t45 * t46 * t16 + t41) * t28 * t26 * t20)
     # - 64 * t5 - 48 * t89 * t27 * (t85 * t1 * t24 * t38 * t12 - t23 * 
     #t18 * t30 * t32 * (t1 * t48 + t23 * t4) + t85 * t4 * t75 * t12) - 
     #128 * t2 - 20 * t88 * t6 * t57 * t73 * t74 * (t52 * t55 + t1) - 12
     # * t10 * t88 * t27 * t56 * t55 * t73 * t74 * (t1 * t58 - t10) - 40
     # * t9

      hjetmass_qqbgg_bubble_pmpm_s12 = -ret/16q0*(0,1q0)
      return

      end function

      complex*32 function hjetmass_qqbgg_bubble_pmpm_s34
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt

      t1 = za(i2, i4)
      t2 = za(i1, i4)
      t3 = zb(i3, i1)
      t4 = zb(i3, i2)
      t5 = t2 * t3
      t6 = t1 * t4
      t7 = t6 + t5
      t8 = zb(i4, i1)
      t9 = za(i1, i3)
      t10 = za(i2, i3)
      t11 = zb(i4, i2)
      t12 = t10 * t11 + t8 * t9
      t9 = t9 * t3
      t13 = t10 * t4
      t2 = t2 * t8
      t11 = t1 * t11
      t14 = t9 + t13 - t2 - t11
      t14 = t14 ** 2
      t15 = 4 * t7 * t12 + t14
      t15 = sqrt(t15)
      t16 = -t15 + t9 + t13 - t2 - t11
      t17 = za(i1, i2)
      t18 = zb(i2, i1)
      t14 = 4 * t7 * t12 + t14
      t14 = sqrt(t14)
      t19 = -t14 + t9 + t13 - t2 - t11
      t20 = 0.1q1 / t12
      t19 = 0.1q1 / t19
      t21 = 2
      t22 = (0.1q1 / 0.2q1)
      t23 = t22 * t16 * t20
      t24 = t21 * t7
      t25 = -t24 * t19 - t23
      t15 = t15 + t9 + t13 - t2 - t11
      t26 = t20 * (t16 - t15)
      t27 = -t23 * t8 + t3
      t28 = 0.1q1 / t12
      t29 = t17 * t18
      t30 = (-t29 - t2 - t11) * t1
      t31 = t10 ** 2 * t4
      t32 = t31 * t16
      t33 = (-t29 + t9) * t28
      t34 = t28 ** 2
      t35 = t28 * t34
      t36 = t16 ** 2
      t37 = t10 * t36 * t34
      t38 = t11 + t2
      t39 = t38 * t28
      t40 = (0.3q1 / 0.4q1) * t37 * t12
      t41 = (t13 + t9) * t28
      t42 = t1 * (t29 + t9)
      t43 = (t29 - t2 - t11) * t28
      t44 = -3 * t13 * t1
      t45 = t1 * t20
      t46 = t45 * (t29 + t2 + t11)
      t38 = t38 * t34
      t47 = t1 * t7
      t48 = t7 * t28
      t49 = -t10 * t16 * t36 * t35 * t12 / 8
      t50 = -t22 * t16 * (t38 * t16 * t10 + t46) - t37 * (t29 - t9 - t13
     #) / 4 + t48 * t16 * t10 + t47 + t49
      t51 = t29 - t2 - t11
      t52 = t45 * (t29 + t9)
      t53 = t9 * t10
      t4 = t1 ** 2 * t4
      t45 = (0.3q1 / 0.2q1) * t45 * t13
      t37 = t22 * t16 * (t53 * t16 * t34 + t32 * t34 + t52) + t37 * t51 
     #/ 4 + t5 * (t16 * t10 * t28 + t1) + t4 + t49 + t45 * t16
      t49 = -t22 * t8 * t15 * t20 + t3
      t54 = t15 ** 2
      t55 = t10 * t54 * t34
      t56 = -t10 * t15 * t54 * t35 * t12 / 8
      t29 = -t22 * t15 * (t38 * t15 * t10 + t46) + t48 * t15 * t10 + t47
     # + t55 * (-t29 + t9 + t13) / 4 + t56
      t2 = t14 + t9 + t13 - t2 - t11
      t2 = 0.1q1 / t2
      t9 = t22 * t15 * t20
      t11 = -t24 * t2 - t9
      t13 = t15 * t10
      t14 = t31 * t15
      t4 = t22 * t15 * (t53 * t15 * t34 + t14 * t34 + t52) + t5 * (t13 *
     # t28 + t1) + t4 + t55 * t51 / 4 + t56 + t45 * t15
      t12 = (0.3q1 / 0.4q1) * t55 * t12
      t20 = 0.1q1 / t25
      t17 = 0.1q1 / t17
      t11 = 0.1q1 / t11
      t22 = 0.1q1 / t26
      t18 = 0.1q1 / t18
      t24 = t29 - t4
      t4 = -t29 + t4
      t25 = t50 - t37
      t26 = t22 ** 2
      t2 = t2 * t11
      t12 = t2 * t15
      t29 = t18 * t35
      t31 = t1 * t8 + t10 * t3
      t35 = t16 * t19 * t20
      t18 = t18 * t34
      ret = 64 * t29 * t26 * t17 * t7 * (t16 * (t27 * (-t50 + t37) * t19
     # * t20 ** 2 + (t25 * t8 + t27 * (-t33 * t16 * t10 + t21 * t10 * (t
     #39 * t16 - t5 - t6) - t32 * t28 - t30 + t21 * t10 * (t41 * t16 + t
     #5) + t43 * t16 * t10 + t42 - t44)) * t19 * t20) + t12 * (t24 * t8 
     #+ t49 * (t4 * t11 - t13 * t33 - t14 * t28 + t21 * t10 * (t39 * t15
     # - t5 - t6) - t30 + t21 * t10 * (t41 * t15 + t5) - t44 + t13 * t43
     # + t42))) + 24 * t18 * t22 * t17 * t7 * (-t35 * (t1 * t27 + t3 * (
     #-t23 * t10 - t1)) + t12 * (t1 * t49 + t3 * (-t9 * t10 - t1))) + 25
     #6 * t29 * t22 * t26 * t17 * t7 * (t12 * t49 * t4 + t35 * t25 * t27
     #) - 128 * t18 * t26 * t17 * t7 * (t25 * t19 * t20 * t27 + t2 * t49
     # * t24) - 4 * t29 * t22 * t17 * t7 * (-t36 * t19 * t20 * t31 + t2 
     #* t31 * t54)

      hjetmass_qqbgg_bubble_pmpm_s34 = -ret/16q0*(0,1q0)
      return

      end function

      complex*32 function hjetmass_qqbgg_bubble_pmpp_mhsq
     &     (i1,i2,i3,i4,za,zb,mt,p,flip)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          double precision mt
          real*16 p(mxpart,4)
          complex*32 hjetmass_qqbgg_bubble_pmpp_s124
          complex*32 hjetmass_qqbgg_bubble_pmpp_s123
          complex*32 hjetmass_qqbgg_bubble_pmpp_s12
          logical flip


      hjetmass_qqbgg_bubble_pmpp_mhsq = 
     & -hjetmass_qqbgg_bubble_pmpp_s124(i1,i2,i3,i4,za,zb,mt,p,flip)
     & -hjetmass_qqbgg_bubble_pmpp_s123(i1,i2,i3,i4,za,zb,mt,p,flip)
     & -hjetmass_qqbgg_bubble_pmpp_s12(i1,i2,i3,i4,za,zb,mt)
      return

      end function

      complex*32 function hjetmass_qqbgg_bubble_pmpp_s123
     &     (i1,i2,i3,i4,za,zb,mt,p,flip)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt
          real*16 p(mxpart,4)
          complex*32 alpha
          logical flip

      alpha = (za(i1,i2)*zb(i2,i1) + za(i1,i3)*zb(i3,i1) +
     & za(i2,i3)*zb(i3,i2))/(za(i1,i2)*zb(i2,i1) + za(i1,i3)*zb(i3,i1))

      p(5,:) = real(alpha)*p(i1,:)
      p(6,:) = (1q0-real(alpha))*p(i1,:) + p(i2,:) + p(i3,:)
      if (flip .eqv. .true.) then
          call spinoru_qp(6,p,zb,za)
      else
          call spinoru_qp(6,p,za,zb)
      end if

      t1 = za(5, i4)
      t2 = za(6, i2)
      t3 = zb(i3, 5)
      t4 = zb(i3, i1)
      t5 = za(6, i4)
      t6 = za(i1, i2)
      t7 = za(i2, i4)
      t8 = zb(i2, i1)
      t9 = za(i2, i3)
      t10 = t6 * t8
      t11 = za(i1, i3) * t4
      t12 = t9 * zb(i3, i2)
      t13 = t11 + t12 + t10
      t14 = za(5, i2)
      t15 = 2
      t16 = t14 * t5
      t17 = t1 * t2
      t18 = -t17 * t15 + t16
      t19 = zb(i1, 6)
      t20 = zb(i2, 5)
      t21 = za(6, i3)
      t22 = zb(i4, 5)
      t23 = t2 * t20
      t24 = t21 * t3
      t25 = t5 * t22
      t26 = t25 + t23 + t24
      t27 = za(6, i1)
      t28 = zb(i2, 6)
      t29 = za(5, i3)
      t30 = zb(i3, 6)
      t31 = zb(i4, 6)
      t32 = t29 * t30
      t33 = t1 * t31
      t34 = t14 * t28
      t35 = t33 + t32 + t34
      t36 = t27 * t19
      t37 = t14 * t20
      t38 = t2 * t28
      t39 = t29 * t3
      t40 = t21 * t30
      t41 = t1 * t22
      t42 = t5 * t31
      t43 = t42 - t41 - t39 + t40 - t37 + t38 + t36
      t43 = t43 ** 2
      t44 = 4 * t26 * t35 + t43
      t44 = sqrt(t44)
      t45 = t42 - t41 + t44 - t39 + t40 - t37 + t38 + t36
      t46 = 0.1q1 / t35
      t47 = (0.1q1 / 0.2q1)
      t48 = t45 * t46
      t49 = t48 * t47
      t50 = t49 * t14 - t2
      t51 = zb(i4, i3)
      t43 = 4 * t26 * t35 + t43
      t43 = sqrt(t43)
      t52 = t42 - t41 + t43 - t39 + t40 - t37 + t38 + t36
      t52 = 0.1q1 / t52
      t53 = t15 * t26
      t54 = t53 * t52 + t49
      t44 = t42 - t41 - t44 - t39 + t40 - t37 + t38 + t36
      t55 = t46 * (-t45 + t44)
      t56 = zb(i4, i1)
      t57 = t7 * zb(i4, i2)
      t58 = t51 * za(i3, i4)
      t59 = za(i1, i4) * t56
      t60 = t2 * t31
      t61 = t14 * t22
      t62 = t2 * t30
      t63 = t3 * t14
      t64 = (t59 + t10 + t11 + t12 + t36 - t37 + t38 + t57 + t58) * t2
      t65 = (t61 + t60) * t5
      t66 = (t63 + t62) * t21
      t67 = -t37 + t38 + t36
      t68 = 4 * t20
      t67 = t34 * t68 * t2 + t67 ** 2
      t67 = sqrt(t67)
      t69 = -t67 - t37 + t38 + t36
      t67 = t67 - t37 + t38 + t36
      t70 = t26 * t2
      t71 = t44 * t46
      t72 = t71 * t47
      t73 = t72 * t14 - t2
      t74 = t72 * t30 + t3
      t43 = t42 - t41 - t43 - t39 + t40 - t37 + t38 + t36
      t43 = 0.1q1 / t43
      t53 = t53 * t43 + t72
      t75 = t49 * t30 + t3
      t76 = 0.1q1 / t1
      t77 = t16 * t76
      t78 = -t2 + t77
      t79 = t5 * t76
      t80 = t79 * t30 + t3
      t81 = 0.1q1 / t31
      t82 = t22 * t81
      t83 = t79 + t82
      t84 = -t5 * t78
      t85 = 0.1q1 / t35
      t86 = t33 + t32
      t87 = t85 ** 2
      t88 = t87 ** 2
      t89 = t85 * t87
      t90 = t46 * ((t59 + t10 + t11 + t12 + t36 - t37 + t38 + t57 + t58)
     # * t2 + (t63 + t62) * t21 + (t61 + t60) * t5)
      t91 = t2 * t86
      t92 = t91 * t87
      t93 = -t59 - t10 - t11 - t12 + t36 - t38 + t40 + t42 - t57 - t58
      t94 = t44 ** 2
      t95 = t44 * t94
      t96 = t87 * t94
      t97 = t96 * t14
      t98 = t41 + t39
      t99 = t98 * t85
      t100 = (0.1q1 / 0.4q1)
      t101 = -t47 * t44 * (-t92 * t44 + t90) - t100 * t97 * t93 + t2 * (
     #t99 * t44 - t23 - t24 - t25)
      t72 = t72 * t31 + t22
      t49 = t49 * t31 + t22
      t102 = t61 * t81 + t2
      t103 = t22 * t30
      t104 = t103 * t81 - t3
      t105 = 0.1q1 / t14
      t106 = 0.1q1 / t28
      t107 = t69 * t105
      t108 = t107 * t47 * t106
      t109 = -t108 - t82
      t110 = t67 * t105
      t111 = t110 * t47 * t106
      t112 = -t82 - t111
      t113 = -t108 + t79
      t114 = t79 - t111
      t115 = t69 * t106
      t116 = t115 * t47 - t2
      t117 = t108 * t31 + t22
      t118 = 0.1q1 / t69
      t119 = t23 * t15
      t120 = t119 * t118 + t108
      t121 = t106 * t105
      t122 = t121 * (t69 - t67)
      t123 = t107 * t106
      t124 = t123 - t71
      t125 = t123 - t48
      t126 = t67 * t106
      t127 = t126 * t47 - t2
      t128 = t111 * t31 + t22
      t129 = 0.1q1 / t67
      t119 = t119 * t129 + t111
      t130 = t110 * t106
      t71 = t130 - t71
      t48 = t130 - t48
      t108 = t108 * t30 + t3
      t111 = t111 * t30 + t3
      t131 = t20 * t106
      t132 = t2 * t105
      t133 = -t132 + t131
      t134 = t10 * t106
      t135 = t134 + t2
      t136 = t17 * t105
      t137 = -t136 + t5
      t138 = t5 * t135
      t139 = t2 * t5
      t140 = t36 * t106
      t141 = t105 ** 2
      t142 = t106 ** 2
      t143 = t106 * t142
      t144 = t67 ** 2
      t145 = t67 * t144
      t146 = t5 * t20
      t147 = t10 * t105
      t148 = t147 * (t136 + t5)
      t149 = t142 * t144
      t150 = t36 - t10
      t151 = t121 * t1
      t152 = t151 * t150 + t5
      t153 = t142 * t105
      t154 = t153 * t145
      t155 = t2 ** 2
      t156 = t155 * t20
      t157 = (0.1q1 / 0.8q1)
      t158 = t100 * t149 * (-t17 * t67 * t141 + t146 + t148) - t47 * t11
     #0 * (t2 * (t126 * t133 * t1 + t138) + t140 * (t126 * t137 + t139))
     # + t157 * t154 * t152 - t156 * (-t130 * t1 + t5)
      t159 = t69 ** 2
      t160 = t69 * t159
      t161 = t17 * t69 * t141
      t162 = t153 * t160
      t137 = t100 * t142 * t159 * (t148 + t146 - t161) + t157 * t162 * t
     #152 - t47 * t107 * (t2 * (t115 * t133 * t1 + t138) + t140 * (t115 
     #* t137 + t139)) - t156 * (-t123 * t1 + t5)
      t138 = t45 ** 2
      t148 = t45 * t138
      t163 = t87 * t138
      t164 = t163 * t14
      t90 = -t100 * t164 * t93 - t47 * t45 * (-t92 * t45 + t90) + t2 * (
     #t99 * t45 - t23 - t24 - t25)
      t92 = (-t59 - t10 - t11 - t12 + t36 - t38 + t40 + t42 - t57 - t58)
     # * t14 * t85
      t93 = t86 * t85
      t99 = t15 * t2 * (t93 * t45 + t39 + t41) - t92 * t45 - t64 - t65 -
     # t66
      t165 = -t136 - t5
      t166 = t2 * (t38 + t36)
      t167 = t36 * t105 - t20
      t168 = t107 + t20
      t169 = t36 * t115
      t170 = t106 * t159
      t171 = (0.3q1 / 0.2q1)
      t172 = (0.3q1 / 0.4q1)
      t173 = t121 * t17
      t174 = t173 * t171
      t175 = t110 + t20
      t176 = t175 * t155
      t177 = t106 * t144
      t178 = t14 ** 2
      t179 = t5 ** 2
      t180 = t76 ** 2
      t181 = t79 * t28
      t182 = t14 * t2
      t183 = -t15 * t182 * t179 * t76 * (t20 + t181) + t5 * (t178 * t5 *
     # t179 * t28 * t76 * t180 + t156 + t79 * t2 * (t38 + t36) + t14 * (
     #t37 - t36) * t180 * t179)
      t179 = t179 * t178
      t184 = t36 * t2
      t185 = t10 * t2
      t186 = t1 * t155
      t187 = -t115 + t2
      t188 = t42 - t10 + t40
      t189 = t20 ** 2
      t190 = t19 ** 2
      t191 = t121 * t86
      t192 = t131 * t22
      t193 = t155 * t141
      t194 = t59 * t106
      t195 = t2 * t189
      t196 = t123 * t2
      t197 = t147 * t142
      t27 = t190 * t27 ** 2
      t198 = t42 - t39 + t40 - t41
      t199 = -t33 - t32
      t200 = t10 * t86
      t201 = t36 * t86
      t202 = -t201 + t200
      t203 = t106 * (t121 * t202 + t11 + t12 - t40 - t42 + t57 + t58 + t
     #59) + t2
      t204 = t142 * t141
      t205 = t36 * t121
      t206 = t10 * t121
      t207 = t143 * t105
      t146 = t207 * t160 * (t20 * (-t59 - t11 + t40 - t57 - t58) + (-t32
     # * t69 + t200) * t141 * t2) + t162 * (t106 * (t31 * (-t161 + t146)
     # - t12 * t20) + t24 + t25 + t205 * (t59 + t69 + t11 + t12 + t39 + 
     #t41 + t57 + t58) + t206 * (t42 - t69 - t39 + t40 - t41))
      t161 = t132 * t199
      t162 = t25 + t24
      t208 = t131 * t31
      t209 = t131 * t30
      t210 = t131 * (t161 + t39 + t41)
      t211 = t139 * (t208 + t22)
      t212 = t2 * (t209 + t3) * t21
      t213 = -t41 + t10 - t39
      t214 = t10 * t162
      t215 = t36 * t162
      t216 = t196 * (-t206 * t159 + t214 + t215 + t107 * (t2 * (-t39 + t
     #10) + t140 * t213)) + t107 * (t2 * (t11 * t131 * t187 + t12 * t131
     # * t187 + t57 * t131 * t187 + t58 * t131 * t187 + t59 * t131 * t18
     #7 + t211 + t212 + t115 * (-t136 * t22 + t210)) + t169 * (t105 * t1
     #55 + t106 * t162))
      t217 = (0.1q1 / 0.16q2)
      t123 = t100 * (t170 * t141 * t2 * (t2 * (t58 + t11 + t40 + t57) + 
     #t134 * t198 + t115 * (t29 * (-t132 * t30 + t3) - t11 - t57 - t58 -
     # t59) + t140 * (t123 * t199 + t11 + t12 + t40 + t42 + t57 + t58 + 
     #t59)) + t159 * (t195 * t142 + t36 * t107 * t143 * (t105 * t188 - t
     #20) + t12 * t2 * t141 * t106 * t187 + t196 * (t151 * t22 + t131 * 
     #(t191 + 1) - t132 * (t151 * t31 + 1)) + t197 * (t115 * t20 - t24 -
     # t25) + t106 * (t193 * t31 - t192) * t5 + t193 * (t194 + t2) + t27
     # * t142 * t141 * (t115 + t2) - t24 * t20 * t142)) - t157 * t146 + 
     #t47 * t216 + t217 * t204 * t159 ** 2 * t203 + t23 * (t2 * (t2 * (t
     #107 + t20) + t24 + t25 + t115 * (t105 * (-t41 - t69 + t36 - t39) -
     # t20)) + t10 * t123 * t187)
      t140 = -t126 + t2
      t146 = t147 + t20
      t147 = t126 - t2
      t151 = t130 * t2
      t187 = t121 * t144
      t196 = t204 * t144 ** 2 * t203
      t216 = -t59 - t57 - t58
      t218 = (t42 - t11 - t12 + t40) * t20
      t207 = t207 * t145 * ((-t67 * t86 + t200) * t141 * t2 + t218) + t1
     #54 * (t131 * t216 + t24 + t25 + t205 * (t59 + t67 + t11 + t12 + t3
     #9 + t41 + t57 + t58) + t206 * (t42 - t67 - t39 + t40 - t41))
      t98 = t151 * (-t206 * t144 + t214 + t215 + t110 * (t36 * (t106 * t
     #213 + t2) + t185)) + t110 * (t2 * (t126 * (-t132 * t98 + t210) + t
     #211 + t212) + t58 * t23 * t106 * t140 + t194 * t23 * t140 + t12 * 
     #t23 * t106 * t140 + t11 * t23 * t106 * t140 + t57 * t23 * t106 * t
     #140 + t215 * t67 * t142)
      t98 = -t100 * (t187 * (t155 * (t105 * (-t42 - t11 - t12) + t33 * t
     #126 * t141) + t2 * (t204 * t201 * t67 + t121 * (t67 * (-t41 + t11 
     #+ t12) - t27 + t36 * (-t59 - t11 - t12 - t40 - t42 - t57 - t58) + 
     #t10 * (-t42 + t39 - t40 + t41))) + t20 * t67 * t142 * t150) + t144
     # * (-t2 * (t142 * t189 + t193) + t36 * t67 * t141 * t143 * (-t42 +
     # t10 - t40) - t27 * t67 * t141 * t143 + t58 * t2 * t141 * t106 * t
     #147 + t194 * t2 * t141 * t147 + t57 * t2 * t141 * t106 * t147 + t1
     #51 * (-t39 * t121 + t131 * (-t191 - 1) + t132 * (t32 * t121 + 1)) 
     #+ t106 * (t146 * t106 * t3 - t193 * t30) * t21 + t25 * t142 * t146
     #)) - t157 * t207 + t217 * t196 + t47 * t98 + t23 * (t2 * (t2 * (t1
     #10 + t20) + t24 + t25 + t126 * (t105 * (-t41 - t67 + t36 - t39) - 
     #t20)) + t10 * t130 * t140)
      t130 = t25 + t23 + t24
      t140 = t38 + t10
      t143 = t130 * t2
      t146 = (-t41 + t10 + t36 - t37 + t38 - t39) * t2 * t85
      t147 = t28 ** 2
      t151 = t32 * t20
      t191 = t2 * t29
      t193 = t14 * t21
      t194 = t27 * t14
      t196 = t155 * t147
      t206 = t27 * t2
      t207 = t59 + t11 + t12 + t40 + t42 + t57 + t58
      t210 = t33 * t20
      t211 = t91 * t28
      t212 = t178 * t20
      t219 = (t16 - t17) * t31
      t220 = t212 * t10
      t221 = t155 * t28
      t222 = t221 * t207
      t223 = t212 * t162
      t224 = t2 * t147
      t225 = (t59 + t11 + t12 - t40 - t42 + t57 + t58) * t28
      t226 = t14 * (t225 + t224) + t200 - t201
      t227 = t40 + t11 + t12
      t228 = -t25 - t24
      t229 = t23 * t14
      t230 = t185 * t34
      t231 = t46 * t2
      t46 = t44 * t155 * (t46 * t20 * (-t59 - t42 - t57 - t58) + (t20 * 
     #t86 + t28 * (t41 - t10 - t36 + t39)) * t87 * t44) + t44 * (t230 * 
     #t94 * t89 + (t229 * (t59 + t11 + t12 - t39 - t41 + t57 + t58) + t3
     #6 * (t14 * t228 + t2 * (t41 - t10 + t39))) * t87 * t44 + t231 * (-
     #t23 * t227 - t24 * (t38 + t10 + t36) + t25 * (-t38 - t10 - t36)))
      t232 = t59 + t11 + t12 + t39 + t41 + t57 + t58
      t233 = t36 * t232 * t14
      t234 = t10 * (t14 * t198 + t91)
      t235 = t89 * t95
      t150 = t28 * t150 * t85
      t46 = t100 * (t96 * (t184 * (t59 + t57) + t222 - t223 + (t2 * (t14
     # * (t28 * (-t59 + t37 + t41) + t210) - t211) + t220 + t36 * (t219 
     #- t212)) * t85 * t44 - t10 * t3 * (t193 + t191)) + t96 * (t2 * (t1
     #78 * t189 + t196) + t206 + (t36 * (t30 * (t193 - t191) - t10 * t14
     #) + t194 + t182 * (t28 * (-t58 - t11 - t12 - t38 + t39 - t57) + t1
     #51)) * t85 * t44 + t10 * (t2 * (t42 + t40 - t41) - t16 * t22) + t1
     #84 * (t58 + t11 + t12 + t40 + t42))) - t157 * (t235 * t178 * (t20 
     #* (-t59 + t40 + t42 - t57 - t58) + t150 * t44 + t24 * t28) + t235 
     #* (t14 * (t14 * (-t20 * (t12 + t11) + t25 * t28) - t211 * t85 * t4
     #4) + t233 + t234)) + t217 * t88 * t94 ** 2 * t14 * t226 - t47 * t4
     #6 + t23 * (-t97 * t140 + t146 * t44 + t143)
      t94 = t199 * t2
      t97 = t162 * t28
      t193 = t45 * t2 * (-t10 * t34 * t138 * t89 + t231 * (t20 * t207 + 
     #t97)) + t45 * (t231 * (t215 + t214) + (t2 * (t20 * (t14 * (-t59 - 
     #t11 - t12 + t39 + t41 - t57 - t58) + t94) + t38 * t213) + t36 * (t
     #14 * t162 + t2 * t213 + t221)) * t87 * t45)
      t213 = t89 * t148
      t214 = t213 * t178
      t88 = t100 * (t163 * (t184 * (t59 + t11 + t12 + t42 + t57) + (t2 *
     # (t14 * (t28 * (-t59 + t37 + t39) + t151) - t211) + t220 + t36 * (
     #t14 * (t40 - t10 - t37) - t32 * t2)) * t85 * t45 + t221 * (t59 + t
     #42 + t57 + t58) + t10 * (t22 * (-t16 - t17) + t139 * t31)) + t163 
     #* (t178 * (t20 * t228 + t195) + t2 * (t38 * t227 + t196 + t36 * (t
     #58 + t36 + t40)) + (t219 * t36 + t194 + t182 * (t28 * (-t58 - t11 
     #- t12 - t38 + t41 - t57) + t210)) * t85 * t45 + t10 * (t21 * (-t63
     # + t62) - t39 * t2))) - t157 * (t214 * (t150 * t45 + t218 + t97) +
     # t213 * (t14 * (t38 * t199 * t85 * t45 + t37 * t216) + t233 + t234
     #)) + t217 * t88 * t138 ** 2 * t14 * t226 + t47 * t193 + t23 * (-t1
     #64 * t140 + t146 * t45 + t143)
      t100 = -t192 + t132 * (t208 + t22)
      t138 = t59 + t11 + t12 + t57 + t58
      t140 = t132 * t29 - t21
      t143 = -t2 - t134
      t136 = t136 - t5
      t146 = -t131 * t3 + t132 * (t209 + t3)
      t150 = (-t38 + t37) * t2
      t151 = t23 * t22
      t157 = t23 * t3
      t164 = t59 + t11 + t12 + t40 + t42 + t57 + t58
      t192 = t25 + t24
      t193 = t10 * t192
      t199 = t2 * (t2 * (t164 * t20 + t192 * t28) + t193 + t36 * t192)
      t208 = -t59 - t11 - t12 + t40 + t42 - t57 - t58
      t91 = t131 * t208 + t24 + t25 + t10 * (t91 * t106 * t141 + t121 * 
     #t198) + t205 * t232
      t133 = t131 * t135 + t132 * (-t106 * (t59 + t11 + t12 + t57 + t58)
     # - t2) + t191 * t121 * (t133 * t30 + t3) + t173 * (t133 * t31 + t2
     #2) + t27 * t153 + t36 * t142 * (t105 * (t161 - t10 + t40 + t42) - 
     #t20)
      t135 = -t38 + t37
      t141 = t34 * t22
      t142 = t34 * t3
      t94 = t14 * (t37 * t10 - t196 + t38 * (-t59 - t11 - t12 + t37 - t5
     #7 - t58)) + t194 + t191 * (t135 * t30 + t142) + t17 * (t135 * t31 
     #+ t141) + t36 * (t14 * t188 - t212 + t94)
      t97 = t178 * (t20 * t208 + t97) + t233 + t234
      t135 = -t42 + t39 - t40 + t41
      t161 = t192 * t14
      t173 = t2 * t155 * t147
      t32 = t2 * (t33 + t32)
      t33 = -t38 - t10 - t36
      t188 = t224 * t89
      t191 = t41 - t10 + t39
      t194 = t156 * (t41 + t39 + t37 - t38 - t36 - t10)
      t191 = (t2 * (t20 * (t14 * (t59 + t11 + t12 - t39 - t41 + t57 + t5
     #8) + t32) + t38 * t191) + t36 * (t191 * t2 - t161 - t221)) * t85
      t209 = t38 + t10
      t210 = t96 * t178 * t28
      t211 = t85 * t44
      t71 = 0.1q1 / t71
      t120 = 0.1q1 / t120
      t54 = 0.1q1 / t54
      t8 = 0.1q1 / t8
      t119 = 0.1q1 / t119
      t124 = 0.1q1 / t124
      t125 = 0.1q1 / t125
      t215 = 0.1q1 / t45
      t53 = 0.1q1 / t53
      t122 = 0.1q1 / t122
      t216 = 0.1q1 / t6
      t55 = 0.1q1 / t55
      t217 = 0.1q1 / t9
      t48 = 0.1q1 / t48
      t218 = 0.1q1 / t44
      t7 = 0.1q1 / t7
      t219 = t118 ** 2
      t220 = t118 * t219
      t224 = t127 ** 2
      t226 = t127 * t224
      t227 = t55 ** 2
      t228 = t55 * t227
      t231 = t116 ** 2
      t232 = t116 * t231
      t233 = t122 ** 2
      t234 = t125 ** 2
      t236 = t50 ** 2
      t237 = t73 ** 2
      t238 = t71 ** 2
      t239 = t124 ** 2
      t240 = t129 ** 2
      t241 = t129 * t240
      t242 = t48 ** 2
      t243 = t48 * t242
      t244 = t234 * t215
      t245 = t74 * t53
      t246 = t245 * t237
      t247 = t246 * t72 * t46 * t238 * t43 * t239
      t248 = t247 * t218
      t249 = t75 * t236
      t250 = t249 * t49
      t251 = t250 * t88 * t54 * t234
      t252 = t240 * t119
      t253 = t252 * t242
      t254 = t239 * t234
      t255 = t254 * t219
      t256 = t255 * t120
      t257 = t256 * t123 * t108 * t117
      t258 = t23 * t233
      t259 = t19 * t217
      t260 = t259 * t121
      t261 = t260 * t7
      t262 = t35 ** 2
      t263 = t120 ** 2
      t264 = t218 ** 2
      t265 = t215 ** 2
      t266 = t54 ** 2
      t267 = t49 * t75
      t268 = t244 * t75 * t49 * t88 * t242 * t52
      t269 = t74 * t72
      t270 = t218 * t53
      t17 = t270 * t239 * t43 * t238 * t237 * (t46 * (t31 * t74 + t72 * 
     #(-t245 + t30)) + t269 * (-t15 * (t191 * t44 + t194) + t171 * t96 *
     # t94 - t172 * t96 * t97 + t47 * t95 * t14 * (t14 * (t225 * t89 + t
     #188) + t202 * t89) + t211 * (t36 * (t2 * (t59 + t11 + t12 + t40 + 
     #t42) - t210) + t173 + t222 - t223 + t10 * (t2 * (t42 - t39 + t40) 
     #+ t210 - t24 * t14)) + t32 * t235 * t34 + t199 + (t206 + t195 * t1
     #78 + t10 * (-t16 - t17) * t22 + t184 * (t57 + t58)) * t85 * t44 - 
     #4 * t211 * t229 * t209 - 3 * t230 * t96))
      t95 = t108 * t117
      t96 = t128 * t119
      t202 = t111 * t128
      t210 = t253 * t238
      t67 = t210 * t224 * (t98 * (t111 * (-t96 + t31) + t128 * t30) + t2
     #02 * (-t15 * (t2 * (t1 * (t100 * t67 + t151) + t20 * (t126 * t138 
     #+ t150) + t29 * (t146 * t67 + t157) - t185 * t175) + t36 * (t2 * (
     #t110 * t143 - t23) + t126 * t136 * t22 + t126 * t140 * t3)) + t171
     # * t144 * t133 - t172 * t177 * t91 + t47 * t121 * t145 * t203 + t2
     # * (t25 * t36 + t24 * (t36 + t10)) + t182 * t126 * t189 + t176 * t
     #11 + t176 * t12 + t176 * t57 + t176 * t58 + t176 * t59 + t221 * t1
     #92 + t110 * (t155 * (t42 + t40 + t38) - t149 * t36) + t10 * (t139 
     #* t22 + t154) + t156 * (t42 + t40) + t32 * t204 * t145 + t126 * (-
     #t37 * t192 - t193 + t132 * (-t10 * t135 + t36 * t164 + t27)) + 4 *
     # t23 * t67 * t143 - 3 * t187 * t185))
      t110 = t204 * t227
      t17 = t110 * t26 * (t236 * (-t268 * t266 + t244 * (t88 * (t49 * t3
     #0 + t31 * t75) + t267 * (t171 * t163 * t94 - t172 * t163 * t97 - t
     #15 * (t191 * t45 + t194) - t47 * t148 * t14 * (t14 * (t208 * t89 *
     # t28 - t188) + t89 * (t201 - t200)) - t2 * (t24 * t33 + t25 * t33 
     #- t23 * (t40 + t12)) + t32 * t213 * t34 - t156 * (-t57 - t58 - t59
     # - t42 - t11) - (-t184 * t164 - t173 + t10 * (t135 * t2 + t161)) *
     # t85 * t45 - (t178 * (t192 * t20 - t195) + t2 * (-t38 * t164 - t27
     #)) * t85 * t45 - t214 * t28 * (t36 - t10) - 4 * t229 * t45 * t85 *
     # t209 - 3 * t230 * t163)) * t52 * t242 * t54) + t17) * t85
      t27 = t259 * t7
      t32 = t27 * t216
      t33 = t32 * t8
      t11 = t33 * (t196 * t178 * t3 * t22 * t156 * t130 * t262 * t219 * 
     #t240 * t264 * t265 + t258 * t121 * (t231 * (-t255 * t117 * t108 * 
     #t123 * t263 + t256 * (t123 * (t108 * t31 + t117 * t30) + t95 * (-t
     #15 * (t2 * (t1 * (t100 * t69 + t151) + t20 * (t115 * t138 + t150) 
     #+ t29 * (t146 * t69 + t157) - t185 * t168) + t36 * (t2 * (t107 * t
     #143 - t23) + t115 * t136 * t22 + t115 * t140 * t3)) + t171 * t159 
     #* t133 - t172 * t170 * t91 + t47 * t121 * t160 * t203 + t69 * (t14
     # * (t195 * t106 - t131 * t162) + t206 * t121 + t132 * (t153 * t86 
     #* t159 + t62 * t21) + t205 * (t2 * t207 - t170) + t134 * (t132 * t
     #198 - t24 - t25)) + t197 * t160 + t199 + t107 * (t57 + t58 + t59 +
     # t42 + t38 + t11 + t12) * t155 + 4 * t23 * t69 * t143 - 3 * t185 *
     # t121 * t159))) + t67) * t87 + t17)
      t12 = -0.1q1 / t113
      t13 = 0.1q1 / t13
      t17 = -0.1q1 / t114
      t21 = t12 ** 2
      t24 = t17 ** 2
      t25 = t108 * t231
      t29 = t25 * t137 * t120
      t38 = t252 * t224
      t40 = t50 * t54
      t42 = t121 * t23
      t47 = t4 * t26 * t13
      t57 = t47 * t85
      t58 = t4 * t22
      t59 = t7 * t8
      t67 = t59 * (t259 * (t57 * t228 * (t270 * t73 * t101 * t72 * t43 -
     # t40 * t49 * t90 * t215 * t52) + t42 * (-t38 * t111 * t158 * t24 +
     # t29 * t219 * t21) * t216 * t180 * t122 * t233) + t58 * t2 * t3 * 
     #t70 * t35 * t262 * t216 * t13 * t264 * t265 * (t215 + t218))
      t69 = 0.1q1 / t5
      t86 = -t249 * t54 * t125 * t48 * t52 + t246 * t71 * t43 * t124
      t89 = t53 ** 2
      t91 = t52 ** 2
      t94 = t43 ** 2
      t97 = t26 ** 2
      t92 = t73 * (t15 * t2 * (t93 * t44 + t39 + t41) - t92 * t44 - t64 
     #- t65 - t66)
      t93 = t73 * t101
      t100 = t93 * t31
      t107 = t93 * t72
      t130 = t31 * t50
      t131 = t50 * t49
      t132 = t215 * t52
      t133 = t132 * t54
      t134 = t121 * t26 * t51
      t135 = t216 * t19
      t136 = t111 * t119
      t25 = t25 * t120 * t118 * t124 * t125 - t136 * t224 * t71 * t129 *
     # t48
      t138 = t108 * t137 * t219 * t21
      t140 = t51 * t85
      t37 = t135 * t59 * t42 * t122 * (t140 * t25 + (t231 * (-t138 * t26
     #3 + (t108 * (t15 * (t1 * (t115 * t167 * t2 + t155 * t168) - t169 *
     # t5) + t172 * t170 * t152 - t5 * (-t37 * t115 + t166) - t10 * (t11
     #5 * t165 + t139) - t174 * t159) + t137 * t30) * t21 * t219 * t120)
     # + t38 * t24 * (t111 * (t15 * (t1 * (t126 * t167 * t2 + t176) - t3
     #6 * t126 * t5) + t172 * t177 * t152 - t5 * (-t37 * t126 + t166) - 
     #t10 * (t126 * t165 + t139) - t174 * t144) + t158 * (-t136 + t30)))
     # * t217 * t180 * t122)
      t37 = t37 + t59 * (t135 * (-t212 * t155 ** 2 * t147 * t3 * t69 * t
     #217 * t219 * t240 + t134 * t85 * t55 * t86) + (t259 * t85 * (t43 *
     # (-t107 * t218 * t89 + t270 * (t72 * (t101 * t14 + t92) + t100)) +
     # t133 * (t90 * (t49 * (-t40 + t14) + t130) + t131 * t99)) * t227 *
     # t26 + t35 * t216 * (-t269 * t237 * t94 * t89 * t264 + t250 * t266
     # * t265 * t91) * t55 * t97 + t264 * t265 * t262 * (t216 * (t70 * (
     #t22 * (t63 - t62) - t60 * t3) + t2 * t3 * t22 * (t15 * t2 * (t41 +
     # t39) - t64 - t65 - t66)) + t259 * t2 * t22 * t70)) * t13 * t4)
      t39 = t49 * t236
      t41 = t237 * t72
      t42 = t26 * t55
      t60 = t23 * (t117 * t231 * t120 * t118 * t124 * t125 - t96 * t224 
     #* t71 * t129 * t48) * t122 + t42 * (t41 * t71 * t43 * t124 * t53 -
     # t39 * t54 * t125 * t48 * t52)
      t62 = t224 * t119
      t63 = t62 * t17 * t129
      t64 = t120 * t231
      t65 = t64 * t118 * t12
      t6 = t6 * t56 * t217
      t66 = t121 * t217
      t9 = t9 * t51
      t70 = t3 * t56
      t96 = t13 * t7
      t115 = t96 * t4
      t126 = t115 * t3
      t135 = t58 * t155 * t3 * t13
      t136 = t35 * t218
      t6 = t8 * (t19 * (t66 * (t140 * t60 * t216 + t23 * t4 * t76 * t7 *
     # t122 * (t65 - t63)) + t57 * t7 * t55 * (t40 * t52 * (t6 - t51) + 
     #(-t6 + t51) * t53 * t43 * t73)) + t262 * (-t135 * t215 * t216 * t2
     #64 * t7 - t135 * t265 * t216 * t218 * t7) + t136 * t2 * (t141 * t2
     #16 * (t58 - t70) * t129 * t217 * t118 * t2 + t126 * (t216 * (t9 + 
     #t61) - t56)) * t215)
      t57 = 0.1q1 / t83
      t61 = 0.1q1 / t113
      t113 = 0.1q1 / t22
      t114 = 0.1q1 / t114
      t135 = -0.1q1 / t109
      t137 = -0.1q1 / t112
      t139 = t108 * t7
      t135 = t139 * t81 * t135 + t217
      t143 = t114 ** 2
      t144 = t61 ** 2
      t145 = t78 ** 2
      t146 = t57 ** 2
      t147 = t18 * t78
      t148 = t78 * t84
      t149 = t69 * t216
      t150 = t80 * t69
      t151 = t150 * t216
      t96 = t96 * t58
      t152 = t2 * t216
      t154 = t111 * t81
      t156 = t80 * t57
      t2 = t8 * (t69 * (t148 * t96 * (-t259 * t76 + t151) * t81 * t146 +
     # t96 * (t259 * (t76 * (t14 * t84 - t147) - t148 * t69) + t149 * (t
     #84 * (-t14 * t80 - t30 * t78) + t147 * t80)) * t81 * t57 + t152 * 
     #(-t126 * t18 * t69 + t34 * (t217 * (-t58 + t70) + t3 * (t70 * t113
     # - t4) * t7) * t129 * t118 * t2)) + t121 * t76 * t19 * (-t58 * t14
     #5 * t57 * t61 * t81 * t114 * t217 * t7 + t23 * (t65 * t135 - t63 *
     # (t154 * t7 * t137 + t217)) * t216 * t122 * t51) + t82 * t32 * t20
     #4 * t143 * t144 * t57 * t69 * t145 * (t76 * (t183 * (t156 - t30) -
     # t80 * (t15 * t20 * (t179 * t76 + t186) + 3 * t5 * (t28 * (t179 * 
     #t180 + t155) + t184) + t5 * (t77 * (-t36 - t10) + t185) + t16 * t2
     # * (-6 * t181 - t68))) + t150 * t183))
      t10 = t111 * t56 - t128 * t4
      t15 = t119 ** 2
      t9 = t9 * t216
      t9 = t8 * (t152 * t217 * (t85 * (t231 * (t4 * t117 ** 2 * t120 * t
     #219 * t124 * t125 - t95 * t56 * t120 * t219 * t124 * t125) + t38 *
     # t128 * t71 * t48 * t10) * t122 * t20 + t195 * t19 * t76 * t7 * (-
     #t232 * t108 * t263 * t220 * t12 + t226 * t111 * t15 * t17 * t241) 
     #* t122 - t136 * t141 * t155 * t19 * t3 * t7 * t118 * t129 * t215) 
     #+ t27 * t4 * t13 * (-t39 * t266 * t215 * t91 + t41 * t94 * t89 * t
     #218) * t55 * t97 + t42 * (t115 * (t132 * t40 * t75 * (t9 - t56) + 
     #t245 * (-t9 + t56) * t218 * t43 * t73) + t66 * t216 * (t237 * (-t2
     #45 * t56 * t72 * t71 * t43 * t124 * t218 + t270 * t4 * t72 ** 2 * 
     #t71 * t43 * t124) + t39 * t133 * t125 * t48 * (-t4 * t49 + t56 * t
     #75))))
      t16 = t69 ** 2
      t18 = t82 * t149
      t20 = t216 * t13 * t4
      t10 = t8 * (t7 * (t134 * t87 * t55 * (t44 * t237 * t71 * t43 * t12
     #4 * t53 - t45 * t236 * t54 * t125 * t48 * t52) * t217 * t190 + t25
     #9 * (t155 * (t152 * t142 * t69 * t118 * t129 - t136 * t58 * t13 * 
     #t215) - t82 * t151 * t153 * t78 * t183 * t57 * t144 * t76 * t143 +
     # t156 * t204 * (t18 * t144 * t76 * t114 * t143 + t18 * t61 * t144 
     #* t76 * t143) * t183 * t145) + t20 * (t156 * t82 * t148 * t1 * t69
     # * t16 + t136 * (t3 * t31 + t103) * t215 * t155)) + t23 * ((t38 * 
     #t17 * (t10 * t217 + t111 * (t154 * t56 * t137 - t4) * t7) + t64 * 
     #(t4 * (t117 * t217 + t139) - t108 * t135 * t56) * t12 * t219) * t2
     #16 * t76 + t140 * t204 * t190 * t217 * t7 * (t64 * t124 * t125 - t
     #62 * t71 * t48)) * t122)
      t18 = t261 * t8 * t85
      t27 = t52 * t265
      t28 = t127 * t129
      t32 = t260 * t216
      t12 = t59 * (t47 * (t43 * (t93 * t269 * t216 * t264 * t89 + (t216 
     #* (t72 * (t101 * (-t14 * t74 - t30 * t73) - t92 * t74) - t100 * t7
     #4) - t107 * t259) * t264 * t53) + t27 * t54 * (t216 * (t90 * (t49 
     #* (-t30 * t50 + t75 * (t40 - t14)) - t130 * t75) - t267 * t50 * t9
     #9) - t131 * t259 * t90)) * t227 + t32 * (t269 * t73 * t237 * t71 *
     # t94 * t124 * t89 * t218 - t267 * t50 * t236 * t266 * t125 * t215 
     #* t48 * t91) * t55 * t97 + t259 * t23 * t216 * t122 * (t23 * t85 *
     # (t95 * t232 * t263 * t220 * t124 * t125 - t202 * t226 * t71 * t15
     # * t241 * t48) + (t138 * t116 * t120 * t106 + t29 * (-t121 * t219 
     #* t12 * t21 - t21 * t220) + t252 * t24 * t158 * t111 * t127 * (t10
     #6 * (-t127 * t105 * t17 + 1) - t28)) * t180 * t122))
      t14 = t22 ** 2
      t15 = t81 ** 2
      t1 = t59 * t69 * (t32 * t78 * t145 * t80 * t14 * t61 * t114 * t146
     # * t15 + (t155 * (-t216 * t30 + t259) + (-t151 * t1 + t259) * t15 
     #* t146 * t145 * t14) * t13 * t4)
      t13 = 0.1q1 / t109
      t14 = 0.1q1 / t112
      t17 = 0.1q1 / t83 ** 2
      t19 = t19 * t51
      t21 = t19 * t76
      t5 = t121 * t216 * t8 * (t5 * t102 ** 2 * t104 * t180 * t7 * t17 *
     # t13 * t14 * (t56 * t104 * t113 - t19 * t81) + (t80 * t15 * t7 * (
     #t150 * t56 - t21) * t114 * t61 * t146 + t81 * (t217 * (t69 * (-t4 
     #* (t79 * t31 + t22) + t56 * t80) - t21) - t150 * t4 * t7) * t114 *
     # t61 * t57) * t145 * t22)
      t7 = t47 * t59 * t227 * t216 * (t107 * t245 * t43 * t264 * (t136 +
     # t55) + t40 * t267 * (t35 * t215 * t265 * t52 - t27 * t55) * t90)
      t13 = t33 * (t26 * (t110 * (t27 * t251 * t242 + t247 * t264) - t15
     #3 * t85 * (t245 * t73 * t72 * t46 * t238 * t43 * t239 * t218 + t40
     # * t268) * t227) + t258 * t87 * (t64 * t95 * t254 * t123 * t220 - 
     #t257 * t116 * t106 + t210 * t202 * t127 * t98 * (t28 - t106)))
      ret = -16384 * t261 * t216 * t8 * t85 * (t121 * (t227 * (t248 * (t
     #71 + t124) + (t125 * t234 * t215 * t242 * t52 + t244 * t243 * t52)
     # * t54 * t88 * t49 * t75 * t236) + t228 * (t251 * t215 * t242 * t5
     #2 - t248)) * t26 + t258 * t85 * (t257 * t231 * (-t124 - t125 - t12
     #2) + (-t253 * t71 * t238 + (t122 * t242 - t243) * t240 * t119 * t2
     #38) * t98 * t111 * t128 * t224)) - 4096 * t11 + 1024 * t67 - 256 *
     # t37 - 64 * t6 + 16 * t2 - 128 * t9 + 32 * t10 + 320 * t18 * t4 * 
     #t60 - 512 * t12 + 192 * t18 * t56 * (t23 * t25 * t122 + t42 * t86)
     # + 8 * t1 - 2048 * t7 + 8192 * t13 - 4 * t5 - 24 * t186 * t20 * t5
     #9 * t3 * t16

      hjetmass_qqbgg_bubble_pmpp_s123 = -ret/16q0*(0,1q0)
      return

      end function

      complex*32 function hjetmass_qqbgg_bubble_pmpp_s124
     &     (i1,i2,i3,i4,za,zb,mt,p,flip)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt
          real*16 p(mxpart,4)
          complex*32 alpha
          logical flip

      alpha = (za(i1,i2)*zb(i2,i1) + za(i1,i4)*zb(i4,i1) +
     & za(i2,i4)*zb(i4,i2))/(za(i1,i2)*zb(i2,i1) + za(i1,i4)*zb(i4,i1))

      p(5,:) = real(alpha)*p(i1,:)
      p(6,:) = (1q0-real(alpha))*p(i1,:) + p(i2,:) + p(i4,:)
      if (flip .eqv. .true.) then
          call spinoru_qp(6,p,zb,za)
      else
          call spinoru_qp(6,p,za,zb)
      end if

      t1 = za(5, i2)
      t2 = zb(i1, 6)
      t3 = za(6, i2)
      t4 = zb(i2, 5)
      t5 = za(6, i3)
      t6 = zb(i3, 5)
      t7 = za(6, i4)
      t8 = zb(i4, 5)
      t9 = t3 * t4
      t10 = t5 * t6
      t11 = t7 * t8
      t12 = t11 + t9 + t10
      t13 = zb(i2, 6)
      t14 = za(5, i3)
      t15 = zb(i3, 6)
      t16 = za(5, i4)
      t17 = zb(i4, 6)
      t18 = t1 * t13
      t19 = t14 * t15
      t20 = t16 * t17
      t21 = t18 + t19 + t20
      t22 = za(6, i1) * t2
      t23 = t1 * t4
      t24 = t3 * t13
      t25 = t14 * t6
      t26 = t5 * t15
      t27 = t16 * t8
      t28 = t7 * t17
      t29 = t22 - t23 + t24 - t25 + t26 - t27 + t28
      t29 = t29 ** 2
      t30 = 4 * t12 * t21 + t29
      t30 = sqrt(t30)
      t31 = t22 - t23 + t24 - t25 + t26 - t27 + t28 - t30
      t32 = 0.1q1 / t21
      t33 = (0.1q1 / 0.2q1)
      t34 = t33 * t1
      t35 = t34 * t31 * t32 - t3
      t36 = t33 * t15
      t37 = t36 * t31 * t32 + t6
      t38 = za(i2, i4)
      t39 = zb(i4, i3)
      t40 = za(i1, i2)
      t41 = zb(i2, i1)
      t42 = zb(i3, i1)
      t43 = za(i2, i3)
      t44 = za(i1, i4)
      t45 = zb(i4, i1)
      t46 = 0.1q1 / t21
      t47 = t1 * t7
      t48 = t3 * t16
      t49 = t48 + t47
      t50 = t46 ** 2
      t51 = t46 * t50
      t52 = t31 ** 2
      t53 = t31 * t52
      t54 = (t25 + t23) * t7
      t55 = t10 * t16
      t56 = t16 * (t28 - t24) * t46
      t57 = t38 * zb(i4, i2)
      t58 = t57 * t49
      t59 = za(i3, i4) * t39
      t60 = t59 * t49
      t61 = t40 * t41
      t62 = t61 * t49
      t63 = za(i1, i3) * t42
      t64 = t63 * t49
      t65 = t43 * zb(i3, i2)
      t66 = t65 * t49
      t67 = t44 * t45
      t68 = t67 * t49
      t52 = t50 * t52
      t69 = t47 * t8
      t70 = t14 * t3
      t71 = t1 * t5
      t72 = t47 * t17
      t73 = t3 * t7
      t74 = t73 * t32
      t75 = -t71 + t70
      t76 = t7 ** 2
      t77 = t6 * t75 * t7
      t78 = t48 * (t10 + t9)
      t79 = t1 * t8
      t80 = t79 * t76
      t81 = t73 * t12
      t82 = t20 * t73
      t83 = (-t80 + t78 + t77) * t46
      t84 = t19 + t18
      t85 = (t26 - t59 - t63 - t65 - t61 - t67 - t57 + t22) * t16
      t86 = t84 * t7
      t87 = t20 + t19
      t88 = t27 + t25
      t89 = t48 * t52
      t90 = 2
      t91 = (0.3q1 / 0.4q1)
      t92 = -(0.3q1 / 0.8q1)
      t93 = (0.1q1 / 0.4q1)
      t94 = (0.1q1 / 0.8q1)
      t95 = t27 * t73 * t90
      t53 = t33 * t31 * (t74 * (-t28 - t59 - t63 - t65 - t61 - t67 - t57
     # - t22 - t24 - t26) + (t16 * (t3 * (t26 - t23 + t24) + t69) + t7 *
     # (t15 * (-t71 + t70) - t72) + t22 * (t48 - t47)) * t50 * t31) + t9
     #3 * t52 * (t1 * (t56 * t31 + t54 + t55) + t58 + t60 + t62 + t64 + 
     #t66 + t68) + t94 * t53 * t1 * (t85 * t51 + t86 * t51) + t92 * t48 
     #* t53 * t51 * t87 - t91 * t89 * t88 + t83 * t31 + t82 * t52 - t81 
     #+ t95 * t31 * t46
      t29 = 4 * t12 * t21 + t29
      t29 = sqrt(t29)
      t96 = t22 - t23 + t24 - t25 + t26 - t27 + t28 - t29
      t97 = 0.1q1 / t16
      t98 = t7 * t97
      t99 = t33 * t31 * t32
      t100 = t99 - t98
      t96 = 0.1q1 / t96
      t101 = t90 * t12
      t102 = t101 * t96 + t99
      t30 = t22 - t23 + t24 - t25 + t26 - t27 + t28 + t30
      t103 = t32 * (-t31 + t30)
      t34 = t34 * t30 * t32 - t3
      t36 = t36 * t30 * t32 + t6
      t104 = t30 ** 2
      t105 = t30 * t104
      t106 = t104 * t50
      t107 = t1 * (t16 * (-t26 + t59 + t63 + t65 + t61 + t67 + t57 - t22
     #) + t7 * (-t19 - t18))
      t108 = t106 * t48
      t51 = -t33 * t30 * (t74 * (t28 + t59 + t63 + t65 + t61 + t67 + t57
     # + t22 + t24 + t26) + (t16 * (t3 * (-t26 + t23 - t24) - t69) + t7 
     #* (t15 * (t71 - t70) + t72) + t22 * (-t48 + t47)) * t50 * t30) - t
     #91 * t108 * t88 + t92 * t48 * t105 * t51 * t87 + t93 * t106 * (t1 
     #* (t56 * t30 + t54 + t55) + t58 + t60 + t62 + t64 + t66 + t68) - t
     #94 * t107 * t51 * t105 + t82 * t104 * t50 + t83 * t30 - t81 + t95 
     #* t30 * t46
      t29 = t22 - t23 + t24 - t25 + t26 - t27 + t28 + t29
      t54 = t33 * t30 * t32
      t56 = t54 - t98
      t29 = 0.1q1 / t29
      t58 = t101 * t29 + t54
      t60 = t47 * t97
      t62 = t60 - t3
      t64 = t98 * t15 + t6
      t66 = t97 ** 2
      t68 = -t90 * t1 * t3 * t76 * t97 * (t98 * t13 + t4) + t7 * (t3 * (
     #t10 + t9) + t98 * (t3 * (t26 + t22 + t24 - t25) - t10 * t1) + (t1 
     #* (-t26 - t22 + t23 + t25) - t70 * t15) * t66 * t76 + t1 * t84 * t
     #97 * t66 * t7 * t76)
      t74 = 0.1q1 / t17
      t81 = t8 * t74
      t82 = t98 + t81
      t83 = 0.1q1 / t14
      t84 = t71 * t83 - t3
      t92 = t83 ** 2
      t94 = t5 * t83
      t95 = t1 * t16
      t101 = t95 * t5 ** 2
      t49 = t5 * (-t101 * t92 + t94 * t49 - t73)
      t104 = -t98 + t94
      t105 = 0.1q1 / t15
      t109 = t6 * t105
      t110 = t109 + t94
      t111 = t73 * t14
      t112 = t61 + t67 + t57
      t113 = -t70 * t90 + t71
      t114 = t1 * t14 * t76
      t115 = t114 * t66
      t70 = t7 * (-t3 * t5 + t98 * (t71 + t70) - t115)
      t71 = -t5 * t84
      t116 = t94 * t17 + t8
      t117 = t12 * t3
      t118 = (t27 + t25) * t46
      t119 = t3 * t17
      t120 = t1 * t6
      t121 = t3 * t15
      t20 = t20 + t19
      t32 = t32 * (t3 * (t24 + t59 + t63 + t65 + t61 + t67 + t57 + t22 -
     # t23) + t5 * (t121 + t120) + t7 * (t119 + t79))
      t122 = t3 * t20 * t50
      t123 = t28 - t59 - t63 - t65 - t61 - t67 - t57 + t22 - t24 + t26
      t124 = t106 * t1
      t125 = -t33 * t30 * (-t122 * t30 + t32) - t93 * t124 * t123 + t3 *
     # (t118 * t30 - t10 - t11 - t9)
      t126 = t7 * (t119 + t79)
      t127 = t5 * (t121 + t120)
      t128 = t3 * (t24 + t59 + t63 + t65 + t61 + t67 + t57 + t22 - t23)
      t129 = t109 + t98
      t79 = t79 * t74 + t3
      t130 = t81 + t94
      t131 = -t109 + t81
      t132 = t120 * t105 + t3
      t133 = t28 - t24
      t134 = t25 + t23
      t135 = t48 + t47
      t136 = t73 * (-t28 - t59 - t63 - t65 - t61 - t67 - t57 - t22 - t24
     # - t26)
      t55 = (t1 * (t134 * t7 + t55) + t67 * t135 + t65 * t135 + t63 * t1
     #35 + t61 * t135 + t59 * t135 + t57 * t135) * t46
      t23 = (t16 * (t3 * (t26 - t23 + t24) + t69) + t7 * (t15 * t75 - t7
     #2) + t22 * (t48 - t47)) * t46
      t47 = t73 * t16
      t69 = -(0.9q1 / 0.4q1)
      t72 = t99 + t81
      t75 = t54 + t81
      t135 = t3 ** 2
      t137 = t1 ** 2 * t76
      t32 = -t33 * t31 * (-t122 * t31 + t32) - t93 * t52 * t1 * t123 + t
     #3 * (t118 * t31 - t10 - t11 - t9)
      t33 = t54 * t17 + t8
      t54 = t99 * t17 + t8
      t93 = t1 * t46
      t99 = t93 * (t28 - t59 - t63 - t65 - t61 - t67 - t57 + t22 - t24 +
     # t26)
      t20 = t20 * t46
      t118 = -t99 * t30 + t90 * t3 * (t20 * t30 + t25 + t27) - t126 - t1
     #27 - t128
      t20 = -t99 * t31 + t90 * t3 * (t20 * t31 + t25 + t27) - t126 - t12
     #7 - t128
      t27 = t8 * t15
      t99 = t27 * t74 - t6
      t122 = 0.1q1 / t56
      t44 = 0.1q1 / t44
      t123 = 0.1q1 / t40
      t138 = 0.1q1 / t30
      t103 = 0.1q1 / t103
      t112 = 0.1q1 / t112
      t139 = 0.1q1 / t100
      t43 = 0.1q1 / t43
      t58 = 0.1q1 / t58
      t41 = 0.1q1 / t41
      t140 = 0.1q1 / t31
      t102 = 0.1q1 / t102
      t141 = 0.1q1 / t38
      t142 = t34 ** 2
      t143 = t139 ** 2
      t144 = t103 ** 2
      t145 = t102 ** 2
      t146 = t122 ** 2
      t147 = t140 ** 2
      t148 = t1 * t37
      t149 = t15 * t35
      t150 = t37 * t20
      t151 = t150 * t35
      t152 = t17 * t32
      t153 = t45 * t41
      t154 = t40 * t44
      t155 = t35 * t53 * t123 * t143
      t156 = t155 * t41 * t66
      t157 = t37 * t35
      t158 = t157 * t32
      t159 = t158 * t54
      t160 = t35 * t139
      t161 = t34 * t36
      t162 = t142 * t36
      t163 = t162 * t123 * t41
      t164 = t34 * t58
      t165 = t15 * t34
      t166 = t161 * t118
      t167 = t161 * t17 * t125
      t168 = t112 * t138
      t169 = t58 * t138 * t29
      t170 = t43 * t144 * t12
      t171 = t138 ** 2
      t172 = t58 ** 2
      t173 = t1 * t36
      t174 = t161 * t125
      t175 = t174 * t2 * t141
      t174 = t174 * t33
      t176 = t35 * t102
      t177 = t102 * t147
      t178 = t170 * t112
      t20 = t178 * (-t154 * t175 * t29 * t171 * t58 + t153 * (t177 * t12
     #3 * t96 * t54 * t32 * (t37 * (t176 - t1) - t149) + (t58 * (t123 * 
     #(t33 * (t125 * (-t165 - t173) - t166) - t167) - t175) + t174 * t12
     #3 * t172) * t171 * t29)) + t170 * (t96 * (t147 * (t102 * (t112 * (
     #t44 * (t54 * (t32 * (-t149 - t148) - t151) - t152 * t35 * t37) - t
     #153 * t37 * t35 * (t20 * t54 + t152) * t123) + t157 * (-t156 + (-t
     #153 - t154) * t112 * t32) * t141 * t2) + t159 * t44 * t112 * t145)
     # + t156 * t2 * t37 * t141 * t46 * (-t160 + t1) * t102 * t140) + t1
     #69 * (t168 * t44 * (t33 * (t125 * (t36 * (t164 - t1) - t165) - t16
     #6) - t167) + (-t163 * t141 * t46 * t66 * t122 * t146 + t161 * (-t1
     #38 * t34 + t93) * t66 * t141 * t41 * t123 * t146) * t51 * t2))
      t93 = 0.1q1 / t72
      t152 = 0.1q1 / t7
      t156 = 0.1q1 / t104
      t104 = 0.1q1 / t104
      t100 = -0.1q1 / t100
      t167 = 0.1q1 / t82
      t56 = -0.1q1 / t56
      t170 = 0.1q1 / t110
      t175 = 0.1q1 / t75
      t179 = 0.1q1 / t5
      t180 = t100 ** 2
      t181 = t56 ** 2
      t182 = t104 ** 2
      t183 = t35 ** 2
      t184 = t156 ** 2
      t185 = t179 ** 2
      t186 = t179 * t185
      t187 = t62 * t50
      t188 = t70 * t92 * t184
      t189 = t109 * t84
      t81 = t81 * t167
      t190 = t81 * t152
      t191 = t45 * t112
      t192 = t102 * t140
      t193 = t192 * t139 * t96 * t183
      t194 = t2 * t6 * t141
      t195 = t21 * t138
      t196 = t168 * t140
      t27 = -t195 * t194 * t3 * t135 * t152 * t123 * t43 * t41 * t140 + 
     #t196 * t43 * (t44 * (t6 * (t2 * t40 * t141 - t17) - t27) + t153 * 
     #(t123 * (-t17 * t6 - t27) + t194)) * t21 * t135 - t189 * t14 * t71
     # * t116 * t186 * t44 * t43 * t170 * t112 + (t43 * (-t191 * t189 * 
     #t14 * t71 * t116 * t186 * t170 + (t189 * t49 * t83 * t182 * t170 *
     # t179 * (-t104 * t84 + t1) * t66 + t190 * t62 * ((t180 * (t1 * t50
     # * t181 - t187 * t56 * t181) - t187 * t181 * t100 * t180) * t68 * 
     #t64 + t188 * (t156 * t62 + t1)) * t97) * t141 * t2) + t97 * (t169 
     #* (t45 * t36 ** 2 * t175 * t74 * t141 - t42 * t36 * t141 - t42 * t
     #33 * t43) * t122 * t142 + t193 * (-t37 ** 2 * t45 * t93 * t74 * t1
     #41 + t37 * t42 * t141 + t42 * t54 * t43)) * t103 * t12) * t41 * t1
     #23
      t72 = -0.1q1 / t72
      t186 = 0.1q1 / t8
      t187 = 0.1q1 / t130
      t130 = 0.1q1 / t130
      t194 = 0.1q1 / t131
      t197 = 0.1q1 / t129
      t75 = -0.1q1 / t75
      t129 = 0.1q1 / t129
      t198 = t39 * t123
      t199 = t45 * t43
      t200 = -t199 - t198
      t201 = t105 ** 2
      t202 = t2 ** 2
      t203 = t170 ** 2
      t82 = 0.1q1 / t82 ** 2
      t110 = 0.1q1 / t110 ** 2
      t204 = t62 ** 2
      t132 = t132 ** 2
      t79 = t79 ** 2
      t205 = t167 ** 2
      t206 = t74 ** 2
      t207 = t84 ** 2
      t208 = t98 * t79
      t209 = t109 * t207
      t210 = t6 * t45
      t211 = t210 * t186
      t212 = t45 * t64 * t152
      t156 = t83 * t156
      t131 = 0.1q1 / t131
      t213 = t208 * t99
      t214 = t97 * t141
      t215 = t214 * t41
      t94 = t215 * (t2 * (t74 * (-t156 * t8 * t42 * t204 * t167 * t43 * 
     #t105 * t129 + t213 * t198 * t72 * t82 * t75 * t46) + t198 * t204 *
     # t64 * t8 * t205 * t206 * t46 * t100 * t56) + t213 * t45 * t83 * t
     #123 * t186 * t187 * t82 - t10 * t202 * t132 * t39 * t92 * t43 * t1
     #10 * t197 * t201 * t131 * t74 + t209 * t74 * (t200 * t83 * t2 + t4
     #5 * (t94 * t15 + t6) * t123 * t179) * t130 * t170 * t104)
      t72 = t94 + t41 * (t141 * (t123 * (-t97 * (t208 * t99 ** 2 * t45 *
     # t186 * t82 * t46 * t72 * t75 + t209 * t42 * t179 * t104 * t170) +
     # t179 * (t211 - t42) * t152 * t135) + t2 * (t5 * t132 * t42 * t92 
     #* t97 * t43 * t110 * t197 * t105 + t97 * (t6 * t42 * t207 * t104 *
     # t203 * t43 * t201 + t208 * t200 * t74 * t82 * t187) * t83) + t202
     # * (t10 * t207 * t39 * t92 * t97 * t104 * t43 * t203 * t201 * t130
     # * t74 + t11 * t79 * t39 * t83 * t66 * t43 * t105 * t187 * t82 * t
     #194 * t206)) + t204 * (t190 * t123 * t42 * (t156 * t141 + (t64 * t
     #141 + (t98 * t17 + t8) * t43) * t46 * t56 * t100) + t141 * t206 * 
     #t8 * (-t45 * t64 ** 2 * t152 * t100 * t56 * t123 * t46 + t156 * (-
     #t202 * t7 * t39 * t43 * t105 * t129 * t66 + (t199 + t198) * t97 * 
     #t2 - t212 * t123)) * t205))
      t75 = t153 + t154
      t79 = t153 * t123
      t82 = t79 + t44
      t94 = t6 ** 2
      t98 = t112 * t179
      t99 = t3 * t152
      t110 = t154 * t112
      t129 = t179 * t135
      t130 = t123 * t152
      t131 = t130 * t41
      t132 = -t153 - t154
      t186 = t152 ** 2
      t187 = t68 * t64
      t190 = t3 * t123
      t194 = t162 * t122
      t197 = t183 * t37
      t199 = t195 * t135
      t200 = t141 * t123
      t4 = (t204 * (((-t42 * t100 * t56 * t46 * t97 + t188 * t186 * t123
     # + t130 * (t97 * (-t15 * t68 - t64 * (t90 * (t3 * (t16 * (t10 - t1
     #1) - t17 * t76) + t76 * (t97 * (t25 + t28) + t8) * t1 + t16 * t135
     # * t4 + t137 * t4 * t97) - t7 * (-t3 * (t57 + t59 + t63 + t65 + t6
     #1 + t67) + t60 * (t26 + t59 + t63 + t65 + t61 + t67 + t57 + t22)) 
     #+ 3 * t7 * (t3 * (t26 + t22 + t24) + t115 * t15 + t137 * t13 * t66
     #) - 4 * t73 * t134 + t3 * t76 * t97 * (-6 * t18 - 5 * t19))) + t18
     #7 * t152) * t50 * t181 * t180) * t141 * t2 - t130 * t156 * t45) * 
     #t74 * t167 + t187 * t130 * t214 * t2 * t180 * t181 * t74 * t50 * t
     #205) - t190 * (t99 * t195 * t42 * t140 + t191 * t113 * t185)) * t4
     #3 * t8
      t7 = t84 * t113
      t13 = t1 * t71
      t16 = t17 * t84 * t71
      t18 = t2 * t84
      t19 = t18 * t71 * t141
      t5 = t6 * (t84 * t123 * t41 * t179 * (t18 * t49 * t83 * t182 * t14
     #1 * t66 + t191 * t179 * t71 * t116) * t105 * t203 + t179 * (t200 *
     # t2 * t207 * (t90 * (t101 * t83 + t111) - 4 * t48 * t5) * t83 * t1
     #82 * t41 * t66 + t98 * (-t19 * t154 + t153 * (t123 * (t116 * (t7 -
     # t13) - t16) - t19))) * t105 * t170) - t3 * t113 * t8 * t185 * t44
     # * t112 + t131 * (t129 * t45 + t214 * t81 * t2 * t204 * t92 * t184
     # * (t70 * t167 + t90 * (t114 * t97 + t48 * t5) - 4 * t111))
      t4 = t43 * t5 + t41 * (t4 + t200 * ((-t197 * t93 * t96 * t139 * t1
     #02 + t194 * t29 * t58 * t175) * t97 * t46 * t74 * t103 * t39 * t12
     # * t2 + t199 * t6 * t152 * t140 * (t211 - t42))) + t98 * t189 * t4
     #3 * t71 * (t132 * t141 * t83 * t2 + t116 * t44 * t179) * t203 + t1
     #09 * t179 * t43 * (t79 * t207 * t104 * t97 + (t207 * t49 * t182 * 
     #t123 * t41 * t179 * t66 + (-t13 * t132 + t7 * t132) * t112 * t83) 
     #* t141 * t2 + t98 * (t116 * (t7 - t13) - t16) * t44) * t170
      t5 = t138 + t140
      t7 = t3 * t21
      t13 = t40 * t42
      t16 = t13 * t141 + t39
      t18 = t6 * (t7 * t5 - t1)
      t19 = t38 * t39
      t7 = t43 * (t7 * t196 * t8 * (t44 * (t13 + t18 + t19) + t153 * (t1
     #23 * (t19 + t18) + t42)) + (-t215 * t42 * t142 * t122 * t29 * t58 
     #+ t176 * t96 * (-t39 * t112 * t75 + (-t40 ** 2 * t44 * t112 - t153
     # * t40 * t112 + t160 * t41 * t97) * t141 * t42) + t164 * t112 * (t
     #153 * t16 + t154 * t16) * t29) * t46 * t103 * t12 * t2)
      t13 = t3 * t6
      t16 = t13 * t8
      t3 = t117 * (t8 * (-t121 + t120) - t119 * t6) + t16 * (t90 * t3 * 
     #t88 - t126 - t127 - t128)
      t6 = t29 ** 2
      t18 = t21 ** 2
      t22 = t112 * t75
      t23 = (t169 * (t110 * (t125 * (t36 * (-t164 + t1) + t165) + t166) 
     #+ t161 * t41 * (t34 * (t69 * t108 * t87 + t90 * (t23 * t30 + t77 +
     # t78 - t80) + t91 * t124 * (t86 + t85) + t55 * t30 + t136 - 3 * t4
     #8 * t30 * t46 * t88 + 4 * t47 * (t17 * t30 * t46 + t8) + (0.3q1 / 
     #0.2q1) * t95 * t106 * t133) * t146 * t123 * t66 + t191 * (-t125 * 
     #t58 + t118))) + (t102 * (t112 * (t154 * t35 * (t15 * t32 + t150) +
     # t153 * (t32 * (t149 + t148) + t151)) + t183 * t41 * (t15 * t53 + 
     #t37 * ((0.3q1 / 0.2q1) * t95 * t52 * t133 - t91 * t107 * t52 + t55
     # * t31 + t90 * (t23 * t31 + t77 + t78 - t80) - 3 * t48 * t31 * t46
     # * t88 + t69 * t89 * t87 + t136 + 4 * t47 * (t17 * t31 * t46 + t8)
     #)) * t66 * t143 * t123) - t157 * t41 * (t155 * t66 + t191 * t32) *
     # t145) * t140 * t96) * t46 * t144 * t12
      t24 = (t197 * t54 * t96 ** 2 * t147 * t145 - t162 * t33 * t6 * t17
     #1 * t172) * t103 * t12 ** 2
      t25 = t147 * t171
      t26 = t46 * t141 * t2
      t1 = t43 * (t25 * t79 * t112 * t18 * t3 + t26 * (t110 * t192 * t37
     # * t32 * t96 * (-t176 + t1) + (-t163 * t51 * t146 * t66 * t172 + (
     #t15 * t142 * t51 * t146 * t123 * t66 + t191 * (t165 + t173) * t125
     #) * t41 * t58) * t138 * t29) * t144 * t12 + t24 * t79 * t21 * t112
     #) + t43 * (t44 * t112 * t21 * (t25 * t21 * t3 + t24) + (t23 + t13 
     #* t147 * t171 * t18 * (t22 * t117 - t190 * t73 * (t11 + t9 + t10) 
     #* t186 * t41)) * t141 * t2)
      t3 = t43 * (t26 * (-t22 * t192 * t158 * t96 + t169 * t161 * (t34 *
     # t51 * t146 * t123 * t41 * t66 + t22 * t125) - t197 * t192 * t53 *
     # t96 * t123 * t143 * t41 * t66) * t103 * t144 * t12 + t25 * t16 * 
     #t112 * t117 * t21 * t18 * (t44 * t5 + t79 * t5))
      t5 = t195 + t103
      t9 = t42 * t75
      t10 = t43 * t103 * t12
      t6 = t10 * (-t162 * t168 * t2 * t12 * t141 * t75 * t172 * t6 + t19
     #2 * t96 * t35 * ((-t197 * t215 * t96 * t123 * t139 * t102 + t22 * 
     #t157 * t141 * t102 * t96) * t12 * t2 + t112 * t54 * (t19 * t82 + t
     #9)) + t215 * t2 * t12 * t34 * t142 * t36 * t123 * t6 * t138 * t122
     # * t172 + t164 * t168 * t33 * (t19 * (-t79 - t44) - t9) * t29)
      t9 = t10 * t198 * t2 * t97 * t41 * t46 * (t96 * t183 * t139 * t102
     # - t142 * t29 * t122 * t58)
      ret = 512 * t20 + 32 * t27 + 8 * t43 * ((t94 * (t104 * t123 * t203
     # * t41 * t179 * t201 * t97 * t84 * t207 - t98 * t75 * t201 * t203 
     #* t207) + t129 * (t41 * (t99 * t123 - t191) - t110) - t131 * (t64 
     #* t100 * t56 * t46 + t156) * t206 * t205 * t62 * t204 * t8 ** 2) *
     # t141 * t2 + t98 * (t82 * t201 * t179 * t203 * t116 * t207 * t14 *
     # t94 + t135 * t17 * t82)) - 16 * t4 - 64 * t7 + 24 * t14 * t135 * 
     #t8 * t185 * t43 * t112 * (t79 + t44) + 256 * t1 - 1024 * t3 + 2048
     # * t178 * (t174 * t58 * t171 * t29 * (t44 * t5 + t79 * t5) + t159 
     #* (t21 * t82 * t102 * t140 * t147 - t177 * t103 * t82) * t96) - 12
     #8 * t6 - 160 * t10 * t79 * t97 * (-t194 * t169 + t193 * t37) - 48 
     #* t9 - 4 * t72 + t41 * t43 * t123 * (-80 * t199 * t210 * t152 * t1
     #40 + t81 * (12 * t39 * t2 * t97 + 20 * t212) * t46 * t56 * t100 * 
     #t204)

      hjetmass_qqbgg_bubble_pmpp_s124 = -ret/16q0*(0,1q0)
      return

      end function

      complex*32 function hjetmass_qqbgg_bubble_pmpp_s12
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt

      t1 = za(i1, i2)
      t2 = za(i2, i3)
      t3 = zb(i3, i1)
      t4 = za(i2, i4)
      t5 = zb(i4, i1)
      t6 = t2 * t3
      t7 = t4 * t5
      t8 = t6 + t7
      t9 = za(i1, i3)
      t10 = zb(i3, i2)
      t11 = za(i1, i4)
      t12 = zb(i4, i2)
      t13 = t11 * t12
      t14 = t9 * t10
      t15 = t13 + t14
      t16 = t9 * t3
      t17 = t2 * t10
      t18 = t11 * t5
      t19 = t4 * t12
      t20 = t16 - t17 + t18 - t19
      t20 = t20 ** 2
      t21 = 4 * t8 * t15 + t20
      t21 = sqrt(t21)
      t22 = t21 - t16 + t17 - t18 + t19
      t23 = 0.1q1 / t15
      t24 = (0.1q1 / 0.2q1)
      t25 = t24 * t12
      t26 = t25 * t22 * t23 + t5
      t27 = za(i3, i4)
      t28 = zb(i4, i3)
      t29 = 0.1q1 / t15
      t30 = t9 * t22
      t31 = t3 * t12
      t32 = t5 * t10
      t33 = t32 + t31
      t34 = t27 * t28
      t35 = t34 * t10
      t36 = t34 * t3
      t37 = t33 * t4
      t33 = (t11 * t33 + t35) * t29
      t38 = 2
      t15 = 4 * t8 * t15 + t20
      t15 = sqrt(t15)
      t20 = t15 - t16 + t17 - t18 + t19
      t39 = 0.1q1 / t9
      t40 = t2 * t39
      t41 = t24 * t22 * t23
      t42 = t41 - t40
      t43 = 0.1q1 / t10
      t44 = t3 * t43
      t45 = t41 + t44
      t20 = 0.1q1 / t20
      t46 = t38 * t8
      t47 = t46 * t20 + t41
      t21 = t21 + t16 - t17 + t18 - t19
      t48 = t23 * (t22 + t21)
      t49 = -t18 - t16 + t17
      t50 = t29 * t22
      t51 = t32 + t31
      t52 = t29 ** 2
      t53 = t29 * t52
      t54 = (t4 * t51 - t36) * t23
      t35 = t11 * t51 + t35
      t51 = t22 ** 2
      t55 = t52 * t51
      t56 = (0.1q1 / 0.4q1)
      t57 = t38 * t3 * t8
      t58 = -t50 * t3 * t49 - t24 * t22 * (-t30 * t3 * t10 * t52 + t54) 
     #+ t56 * t55 * t35 - t57
      t59 = zb(i2, i1)
      t15 = t15 + t16 - t17 + t18 - t19
      t60 = 0.1q1 / t11
      t61 = t4 * t60
      t62 = t24 * t21 * t23
      t63 = -t62 - t61
      t64 = 0.1q1 / t12
      t65 = t5 * t64
      t66 = -t62 + t65
      t15 = 0.1q1 / t15
      t46 = -t46 * t15 - t62
      t67 = t41 - t61
      t68 = t41 + t65
      t25 = -t25 * t21 * t23 + t5
      t69 = t9 * t21
      t70 = -t62 - t40
      t71 = -t62 + t44
      t62 = -t62 * t10 + t3
      t41 = t41 * t10 + t3
      t72 = -t40 + t61
      t73 = t40 + t44
      t74 = t40 + t65
      t75 = t61 + t44
      t76 = t61 + t65
      t77 = t4 * t10
      t78 = t11 * t3
      t79 = -t77 + t78
      t80 = t77 * t60 + t3
      t81 = t4 * t80
      t82 = t29 * t21
      t83 = t21 ** 2
      t84 = t52 * t83
      t49 = t82 * t3 * t49 + t24 * t21 * (t69 * t3 * t10 * t52 + t54) + 
     #t56 * t84 * t35 - t57
      t54 = t34 + t16 + t18
      t57 = t34 + t16 + t18
      t85 = t22 * t23
      t86 = t24 * t85 * t57 - t8
      t87 = t34 + t17 + t19
      t85 = t24 * t85 * t87 + t8
      t88 = t34 + t17 + t19
      t89 = t21 * t23
      t57 = -t24 * t89 * t57 - t8
      t87 = -t24 * t89 * t87 + t8
      t89 = t40 * t12 + t5
      t90 = t31 * t43 - t5
      t91 = t32 * t64 - t3
      t9 = t77 * t9
      t77 = t11 * (t34 + t17) + t9
      t92 = t78 * t2
      t93 = (-t16 + t34) * t4
      t23 = (t4 * (-t34 + t16) + t92) * t23
      t94 = t13 * t4
      t95 = t19 + t17 - t18
      t9 = t11 * (t34 + t17) + t9
      t96 = t38 * t4 * t8
      t97 = -t24 * t21 * (-t94 * t21 * t52 + t23) + t56 * t84 * t9 + t82
     # * t4 * t95 - t96
      t9 = -t50 * t4 * t95 + t24 * t22 * (t94 * t22 * t52 + t23) + t56 *
     # t55 * t9 - t96
      t23 = t4 ** 2
      t24 = t60 ** 2
      t56 = t60 * t24
      t6 = t4 * (-t14 * t23 * t24 + t6 + t61 * (t17 - t16))
      t94 = t3 * t5
      t95 = t3 ** 2
      t96 = t43 ** 2
      t98 = t43 * t96
      t99 = 0.1q1 / t63
      t100 = 0.1q1 / t2
      t48 = 0.1q1 / t48
      t27 = 0.1q1 / t27
      t46 = 0.1q1 / t46
      t101 = 0.1q1 / t67
      t47 = 0.1q1 / t47
      t102 = 0.1q1 / t4
      t103 = t54 + t88
      t104 = t48 ** 2
      t105 = t47 ** 2
      t106 = t101 ** 2
      t107 = t99 ** 2
      t108 = t62 * t46
      t109 = t12 * t100
      t110 = t25 * t100
      t111 = t41 * t47
      t112 = t46 * t15
      t113 = t112 * t21
      t114 = t53 * t104 * t8 * t1
      t115 = t46 ** 2
      t116 = t62 * t102
      t117 = t47 * t20
      t54 = t114 * t27 * (t21 * (-t87 * (t116 + t110) * t15 * t115 + t11
     #2 * (t102 * (t10 * t87 + t54 * t62) + t109 * t87)) + t117 * t86 * 
     #t22 * (t102 * (-t111 + t10) + t109)) + t114 * (t113 * (t102 * (t27
     # * (t57 * (-t108 + t10) + t62 * t88) + t82 * t1 * t62 * t97 * t99 
     #* t107 * t100 * t24) + t109 * t57 * t27 + t110 * t27 * (-t46 * t57
     # + t54 + t88)) + (t27 * (t100 * (t103 * t26 + t12 * t85) + t102 * 
     #(t10 * t85 + t103 * t41)) * t47 + t27 * (-t41 * t85 * t102 + (-t86
     # - t85) * t100 * t26) * t105) * t20 * t22 - t111 * t1 * t51 * t9 *
     # t24 * t100 * t102 * t29 * t20 * t101 * t106)
      t85 = t86 + t85
      t86 = t41 * t102
      t85 = t85 * t100 * t26 + t86 * t85
      t57 = -t116 * (t57 + t87) + t110 * (-t57 - t87)
      t87 = 0.1q1 / t45
      t88 = 0.1q1 / t42
      t103 = 0.1q1 / t70
      t109 = 0.1q1 / t71
      t114 = t41 * t9
      t118 = t114 * t24 * t106
      t119 = t26 * t58 * t39 * t87 * t88 * t43
      t120 = t62 * t97 * t107 * t24
      t121 = t25 * t49
      t122 = t121 * t39 * t103 * t43 * t109
      t123 = 0.1q1 / t76
      t63 = -0.1q1 / t63
      t124 = 0.1q1 / t75
      t125 = 0.1q1 / t72
      t67 = -0.1q1 / t67
      t126 = t100 * t5
      t127 = t102 * t3 + t126
      t128 = t63 ** 2
      t129 = t67 ** 2
      t130 = t121 * t15 * t109
      t131 = t117 * t88
      t132 = t131 * t87 * t51
      t133 = t117 * t22
      t134 = t102 * t52 ** 2 * t100
      t135 = t81 * t39 * t124 * t125 * t43
      t136 = t80 * t6 * t52 * t128 * t129
      t137 = t65 * t56 * t123 * t100
      t138 = t137 * t1 * (t136 * (t61 * (t63 + t67) - 1) - t135)
      t9 = t1 * (t8 * (t52 * t27 * (t113 * t127 - t133 * t127) * t48 + t
     #134 * t1 * (t24 * ((t47 * (t9 * t10 + t41 * (-t38 * t4 * (t12 * (-
     #t50 * t11 + t4) + t17 - t18) + t50 * t77 + t92 - t93)) - t114 * t1
     #05) * t106 * t20 * t51 + t112 * t107 * t83 * (-t62 * (-t38 * t4 * 
     #(t12 * (t82 * t11 + t4) + t17 - t18) - t82 * t77 + t92 - t93) + (t
     #108 - t10) * t97)) + (t83 * (t46 * (-t130 * t103 ** 2 + (t109 * (t
     #12 * t49 + t25 * (t38 * t3 * (t10 * (-t69 * t29 - t2) + t16 + t18)
     # - t33 * t21 + t36 - t37)) - t121 * t109 ** 2) * t15 * t103) - t13
     #0 * t103 * t115) + t132 * (-t26 * (t38 * t3 * (t10 * (t30 * t29 - 
     #t2) + t16 + t18) + t33 * t22 + t36 - t37) + (t26 * (t87 + t88 + t4
     #7) - t12) * t58)) * t43 * t39) * t104) + t138)
      t17 = 0.1q1 / t68
      t30 = 0.1q1 / t73
      t33 = 0.1q1 / t66
      t37 = 0.1q1 / t59
      t46 = 0.1q1 / t72
      t47 = t30 ** 2
      t49 = t123 ** 2
      t58 = t64 ** 2
      t69 = t5 ** 2
      t72 = t5 * t69
      t77 = t39 ** 2
      t93 = t39 * t77
      t97 = t24 * t4
      t86 = t1 * ((t133 * (-t44 * t26 ** 2 * t100 * t39 * t87 * t88 + t8
     #6 * t60 * t101 * (t65 * t41 * t17 - t3) + t126 * t26 * t39 * t88) 
     #+ t113 * (t116 * t3 * t99 * t60 + t110 * t39 * t103 * (t44 * t25 *
     # t109 - t5) - t65 * t62 ** 2 * t99 * t33 * t60 * t102)) * t37 * t5
     #2 * t48 * t8 - t1 * t39 * t60 * (t2 * t95 * t77 * t102 * t46 * t47
     # * t96 + t97 * t69 * t100 * t125 * t49 * t58))
      t71 = -0.1q1 / t71
      t45 = -0.1q1 / t45
      t70 = -0.1q1 / t70
      t42 = -0.1q1 / t42
      t75 = 0.1q1 / t75
      t73 = 0.1q1 / t73 ** 2
      t81 = t81 * t60 * t43
      t106 = t1 * t100
      t75 = t60 * t75
      t107 = t73 * t3
      t108 = t2 ** 2 * t28
      t88 = t39 * t88
      t110 = t100 * t60
      t114 = t39 * t103
      t116 = t112 * t83
      t62 = t62 * t99
      t121 = t44 * t102
      t11 = t121 * t77 * (t1 * (t107 * (t90 * t3 * (-t13 * t95 * t96 + t
     #7 + t44 * (-t19 + t18)) * t71 ** 2 * t45 ** 2 * t52 + t75 * (t75 *
     # t3 * (t78 * t43 + t4) - t3)) * t96 + t40 * t47 * (-t89 * (t38 * t
     #94 * (t40 * t11 - t4) + t40 * (t4 * (-t32 - t31) + t36 + t40 * t35
     #)) * t52 * t70 ** 2 * t42 ** 2 + (-t38 * t4 * t3 + t40 * t79) * t2
     #4 * t46 ** 2) * t43) - t108 * t39 * t30 * t29 * t70 * t42)
      t13 = t137 * t23 * t28
      t6 = t1 * (t106 * t102 * t53 * (-t51 * t20 ** 2 * t105 * (t41 * t6
     #0 * t101 + t88 * t26) + (t25 * t39 * t103 + t62 * t60) * t15 ** 2 
     #* t115 * t83) * t48 * t8 ** 2 + t7 * (t106 * t56 * (t39 * (t124 * 
     #(t81 * t125 ** 2 + (t60 * t79 - t3) * t43 * t125) + t81 * t125 * t
     #124 ** 2) + t52 * t128 * t60 * t129 * (-t10 * t6 + t80 * (t38 * (t
     #4 * (t14 * t60 + t12) + t18) * t4 + t4 * (t16 + t34) + t92))) * t6
     #4 * t123 + t106 * t24 ** 2 * (t136 + t135) * t64 * t49) + t53 * t2
     #8 * (t117 * (-t110 * t101 - t88 * t102) * t51 - t116 * (t114 * t10
     #2 + t110 * t99)) * t48 * t8 + t11 - t13 * t29 * t63 * t67)
      t11 = t1 ** 2
      t14 = t116 * t103
      t15 = t29 * t11
      t16 = 0.1q1 / t74
      t18 = -0.1q1 / t66
      t19 = -0.1q1 / t68
      t31 = 0.1q1 / t74
      t32 = 0.1q1 / t76 ** 2
      t34 = t107 * t75
      t35 = t28 * t64
      t36 = t106 * t59
      t41 = t40 * t28
      t66 = t24 * t58
      t65 = t65 * t24 * t123
      t25 = t14 * t25 * t109 + t132 * t26
      t26 = t126 * t1
      t68 = t95 * t90
      t10 = t95 * t90 ** 2 * t77 * t37 * t73 * t96 * t29 * t71 * t45 + t
     #77 * t3 * (t3 * (t29 * t70 * t37 * t42 * t89 ** 2 + t60 * t37 * t4
     #6 * t89) - t41 * t60 * t46) * t96 * t47 - t66 * t69 * t91 ** 2 * t
     #37 * t32 * t29 * t18 * t19 + t13 * t39 * t125 * t124 * t43 - t94 *
     # t77 * t37 * t60 * t46 * (t2 * (t40 * t10 + t3) * t16 * t64 * t102
     # + 1) * t43 * t30 - t106 * t72 * t91 * t24 * t32 * t64 * t58 * t29
     # * t18 * t19 + t53 * (-t28 * t100 * t39 * t43 * t25 + (-t116 * t62
     # * t33 * (t26 + t28) + t111 * (-t26 - t28) * t101 * t20 * t17 * t5
     #1) * t102 * t64 * t60) * t48 * t8 - t97 * t94 * (t61 * t12 + t5) *
     # t39 * t100 * t125 * t37 * t124 * t43 * t123 * t64 + t65 * t37 * (
     #t5 * (-t123 * t29 * t63 * t64 * t67 * t80 ** 2 + t39 * t123 * t64 
     #* t125 * t80) - t3 * t39 * t125) - t68 * t28 * t77 * t73 * t98 * t
     #29 * t71 * t45
      t3 = t1 * t10 + t1 * (t3 * (t77 * (-t1 * t59 * t95 * t28 * t102 * 
     #t73 * t96 ** 2 * t29 * t71 * t45 + t121 * t2 * t60 * t46 * t37 * t
     #30 - t34 * t90 * t37 * t96 + t34 * t28 * t98) + t110 * t1 * t8 * t
     #102 * t53 * t48 * (t117 * t51 * t101 + t116 * t99) - t108 * t93 * 
     #t60 * t102 * t46 * t30 * t43 * t16 * t64) + t4 * (t69 * t39 * t24 
     #* t100 * t125 * t37 * t123 * t64 + t5 * t28 * t39 * t56 * t125 * t
     #49 * t58) + t66 * (-t35 * t91 * t29 * t18 * t19 + (t37 * t91 + t35
     #) * t31 * t39) * t32 * t69 - t36 * t72 * t28 * t24 * t32 * t58 ** 
     #2 * t29 * t18 * t19 + t44 * t29 * t30 * t77 * (-t108 * t1 * t59 * 
     #t77 * t30 * t43 * t102 - t41 * t89 * t30 * t43 - t5 * t89 * t37) *
     # t42 * t70 + t65 * t29 * (-t35 * t36 * t23 * t123 * t24 + t61 * (-
     #t35 * t80 * t123 + t106 * t3) + t3 * t80 * t37) * t67 * t63)
      t4 = t60 * t64
      t4 = t134 * t48 * t59 * t28 * t8 * t11 * (-t117 * t22 * t51 * (t4 
     #* t101 * t17 + t88 * t43 * t87) + t112 * (t114 * t43 * t109 + t4 *
     # t99 * t33) * t21 * t83)
      ret = -t38 * t4 - 64 * t54 + 256 * t53 * t27 * t48 * t104 * t8 * t
     #1 * (t117 * t85 * t22 + t113 * t57) - 128 * t52 * t104 * t8 * t1 *
     # (t117 * (t1 * (t55 * t100 * t48 * t102 * (t119 - t118) + t50 * t1
     #00 * t102 * (-t119 + t118)) + t27 * t85) + t112 * (t1 * (t84 * t10
     #0 * t48 * t102 * (t122 - t120) + t82 * t100 * t102 * (-t122 + t120
     #)) + t27 * t57)) - 32 * t9 + 8 * t86 - 16 * t6 - 12 * t15 * (t2 * 
     #(t89 * t93 * t47 * t70 * t42 * t96 * t102 * t95 + t94 * t93 * t30 
     #* t70 * t42 * t43 * t102) + t126 * (-t7 * t80 * t56 * t49 * t58 * 
     #t63 * t67 + t8 * t39 * t102 * t52 * t48 * (t131 * t51 + t14))) - 2
     #0 * t121 * t15 * t39 * (t68 * t39 * t73 * t96 * t71 * t45 + t8 * t
     #100 * t52 * t48 * t25) - 4 * t3

      hjetmass_qqbgg_bubble_pmpp_s12 = -ret/16q0*(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_triangle_pmmp_s123_mhsq_0
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt
          real*16 cg

      t1 = za(i1, i2)
      t2 = zb(i2, i1)
      t3 = za(i1, i3)
      t4 = zb(i3, i1)
      t5 = za(i2, i3)
      t6 = zb(i3, i2)
      t7 = t1 * t2
      t8 = t3 * t4
      t9 = t5 * t6 + t7 + t8
      t10 = za(i1, i4)
      t11 = zb(i4, i1)
      t12 = za(i2, i4)
      t13 = zb(i4, i2)
      t14 = za(i3, i4)
      t15 = zb(i4, i3)
      t16 = t10 * t11
      t17 = t12 * t13
      t18 = t14 * t15
      t19 = t16 + t17 + t18
      if ( real(t19) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t18 = cg * sqrt(t19 ** 2) + t16 + t17 + t18
      t19 = 0.1q1 / t18
      t8 = -2 * t16 * t9 * t19 + t7 + t8
      t20 = t1 * t13 + t15 * t3
      t21 = -2 * t10 * t9
      t22 = t21 * t15 * t19 - t1 * t6
      t21 = t21 * t13 * t19 + t3 * t6
      t23 = t10 * t8
      t24 = t11 * (t12 * t21 + t23)
      t3 = t3 * t11
      t25 = t5 * t13
      t26 = t25 + t3
      t27 = t1 * t11
      t28 = t5 * t15
      t29 = t27 - t28
      t30 = -2 * t12 * t9 * t11 * t19 + t4 * t5
      t31 = t11 * t8
      t32 = t10 * (t13 * t30 + t31)
      t33 = -2 * t14 * t9 * t11 * t19 - t2 * t5
      t1 = 0.1q1 / t1
      t22 = 0.1q1 / t22
      t4 = 0.1q1 / t4
      t18 = 0.1q1 / t18
      t34 = 0.1q1 / t2
      t24 = 0.1q1 / t24
      t35 = t11 ** 2
      t36 = t11 * t35
      t37 = t12 ** 2
      t38 = t24 ** 2
      t39 = t24 * t38
      t23 = t23 * t36 * t38
      t40 = t35 * t24
      t41 = t8 * t15 * t22
      t42 = t13 * t8
      t43 = t42 * t37 * t35
      t44 = t34 * t22
      t45 = t1 * t5
      t46 = t45 * t18
      t47 = 0.1q1 / t12
      t48 = 0.1q1 / t14
      t33 = 0.1q1 / t33
      t32 = 0.1q1 / t32
      t49 = 0.1q1 / t15
      t50 = t11 * t21
      t51 = t42 - t50
      t32 = t10 * t32
      t33 = t33 * t49 - t32
      t52 = t11 * t24
      t53 = t8 ** 2
      t54 = t22 ** 2
      t55 = mt ** 2
      t56 = t15 ** 2 * t53 * t54
      t57 = t8 * (t22 * t48 - t52)
      t58 = t11 * t55
      t59 = t45 * t34
      t60 = t14 * t35
      t61 = t60 * t38
      t54 = t15 * t54
      t62 = t34 * t19 * t4 * t1
      t6 = 0.1q1 / t6
      t63 = 0.1q1 / t9
      t3 = t25 + t3
      t15 = t45 * t15
      t25 = t32 * t26 * t30
      t32 = t14 * t47
      t64 = t9 ** 2
      t65 = t2 * t4
      t66 = t34 * t18 ** 2 * t4 * t1 * t20
      t19 = t11 * ((-t6 * (t65 + t32) + t45 * (-t32 * t34 - t4)) * t63 *
     # t19 * t29 * t55 + t66 * t8 * t12 * t64 * (t61 * (t16 * (t42 * t12
     # * t24 - 1) - t17) + t54))
      t67 = t5 * t35
      t68 = t35 * t34 * t4
      t1 = t45 * (t22 * (-t11 * t34 - t48 * t5) + t52 * (-t12 * t11 * t4
     # + t5)) * t20 + t68 * (-t11 * t49 + t45) - t17 * t67 * t53 * t4 * 
     #t38 + t31 * t4 * (t12 * (t52 * t1 * t34 * t26 + t67 * t21 * t38) +
     # t44 * (t15 - t11)) + t68 * t25 * t1
      t17 = t36 * t14
      t10 = t64 * (t66 * (t14 * (t10 ** 2 * t11 * t35 ** 2 * t53 * t39 +
     # t24 * t36) - t22 * (t56 + t35)) * t12 + t17 * t66 * t12 * t37 * t
     #53 * t13 ** 2 * t39) + t9 * (t17 * t8 * t38 * (-t16 * t8 * t24 + 1
     #) * t18 * t4 * t20 * t12 - t17 * t18 * t37 * t53 * t13 * t20 * t4 
     #* t39) - t58 * t5 * t47 * t63 * (t59 + t6)
      ret = 24 * t46 * t20 * t9 * (t43 * t4 * t38 + t44 * (t41 - t11) + 
     #(-t40 + t23) * t4 * t12) + 8 * t59 * t4 * ((t12 * (t23 * t51 + t40
     # * (-t42 + t50)) + t35 + t56 + t43 * t51 * t38) * t18 * t9 + t58 *
     # (t11 * t33 * t47 * t30 + t57)) - 32 * t62 * t53 * t12 * t55 * (t5
     #4 * t26 * t48 + t61 * (t11 * t20 - t13 * t29)) + 16 * t11 * (t63 *
     # (-t3 * t6 + t4 * (t5 * (t15 - t11) + (t28 - t27) * t6 * t2) - t59
     # * t3) + t62 * (t11 * (t32 * t29 * t30 * t33 + t29 * t49 - t25) + 
     #t57 * t26 * t12) * t55 - t7 * t60 * t12 * t53 * t20 * t4 * t39 - t
     #41 * t46 * t9 * t34 * t4) - 128 * t19 + 4 * t1 - 64 * t10 + 48 * t
     #62 * t58 * t8 * t29 * (-t52 * t14 + t22) - 12 * t65 * t67 * t12 * 
     #t8 * t20 * t38

      hjetmass_qqbgg_triangle_pmmp_s123_mhsq_0 = ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     &  hjetmass_qqbgg_triangle_pmmp_s123_mhsq_0_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          real*16 cg

      t1 = za(i2, i3)
      t2 = zb(i4, i1)
      t3 = za(i1, i2)
      t4 = za(i2, i4)
      t5 = zb(i2, i1)
      t6 = za(i1, i3)
      t7 = zb(i3, i1)
      t8 = zb(i3, i2)
      t9 = t3 * t5
      t10 = t6 * t7
      t11 = t1 * t8 + t10 + t9
      t12 = za(i1, i4)
      t13 = zb(i4, i2)
      t14 = za(i3, i4)
      t15 = zb(i4, i3)
      t16 = t12 * t2
      t17 = t4 * t13
      t18 = t14 * t15
      t19 = t16 + t17 + t18
      if ( real(t19) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t17 = cg * sqrt(t19 ** 2) + t16 + t17 + t18
      t17 = 0.1q1 / t17
      t18 = -0.2q1 * t4 * t11 * t2 * t17 + t1 * t7
      t9 = -0.2q1 * t16 * t11 * t17 + t10 + t9
      t10 = t2 * t9
      t16 = t12 * (t13 * t18 + t10)
      t19 = t1 * t13 + t2 * t6
      t20 = -0.2q1 * t12 * t11
      t21 = t2 * (t12 * t9 + t4 * (t20 * t13 * t17 + t6 * t8))
      t22 = -0.2q1 * t14 * t11 * t2 * t17 - t1 * t5
      t23 = -t1 * t15 + t2 * t3
      t20 = t20 * t15 * t17 - t3 * t8
      t11 = 0.1q1 / t11
      t24 = 0.1q1 / t4
      t8 = 0.1q1 / t8
      t25 = 0.1q1 / t5
      t26 = 0.1q1 / t3
      t27 = t1 * t26
      t28 = t27 * t25
      t29 = 0.1q1 / t15
      t16 = 0.1q1 / t16
      t30 = 0.1q1 / t14
      t20 = 0.1q1 / t20
      t21 = 0.1q1 / t21
      t7 = 0.1q1 / t7
      t22 = 0.1q1 / t22
      t12 = t12 * t16
      t16 = t2 * t21
      t31 = t30 * t20
      t32 = t14 * t24
      t33 = t17 * t26 * t7 * t25
      ret = 0.64q2 * t1 * t2 * t24 * t11 * (t28 + t8) + 0.8q1 * t28 * t7
     # * t2 * (t9 * (t31 - t16) + t2 * (t22 * t29 - t12) * t24 * t18) + 
     #0.48q2 * t10 * t23 * t17 * t26 * t25 * t7 * (-t16 * t14 + t20) - 0
     #.16q2 * t33 * t2 * (t2 * (-t23 * t29 * (t32 * t18 * t22 + 0.1q1) +
     # t12 * t18 * (t32 * t23 + t19)) + t19 * (-t31 + t16) * t9 * t4) + 
     #0.128q3 * t17 * t11 * t23 * t2 * (t8 * (t5 * t7 + t32) + t27 * (t3
     #2 * t25 + t7)) - 0.32q2 * t33 * t9 ** 2 * t4 * (t19 * t15 * t30 * 
     #t20 ** 2 + t14 * t2 ** 2 * t21 ** 2 * (-t13 * t23 + t2 * (t13 * t3
     # + t15 * t6)))


      hjetmass_qqbgg_triangle_pmmp_s123_mhsq_0_rat = ret/32q0*(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_triangle_pmmp_s124_mhsq_0
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt
          real*16 cg

      t1 = za(i2, i3)
      t2 = zb(i3, i1)
      t3 = za(i1, i3)
      t4 = zb(i4, i1)
      t5 = zb(i4, i2)
      t6 = t1 * t5
      t7 = t3 * t4 + t6
      t8 = za(i1, i2)
      t9 = zb(i2, i1)
      t10 = za(i1, i4)
      t11 = za(i2, i4)
      t12 = t8 * t9
      t13 = t10 * t4
      t14 = t11 * t5 + t12 + t13
      t15 = zb(i3, i2)
      t16 = za(i3, i4)
      t17 = zb(i4, i3)
      t18 = t3 * t2
      t19 = t1 * t15
      t20 = t16 * t17
      t21 = t18 + t19 + t20
      if ( real(t21) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t20 = cg * sqrt(t21 ** 2) + t18 + t19 + t20
      t21 = 0.1q1 / t20
      t12 = -2 * t14 * t18 * t21 + t12 + t13
      t13 = 2 * t3
      t18 = t2 * (t1 * (-t13 * t14 * t15 * t21 + t10 * t5) + t12 * t3)
      t22 = t1 * t9 - t16 * t4
      t13 = t13 * t14 * t17 * t21 - t5 * t8
      t23 = -2 * t1 * t14 * t2 * t21 + t11 * t4
      t24 = t15 * t23
      t25 = t3 * (t12 * t2 + t24)
      t26 = 2 * t14 * t16 * t2 * t21 - t11 * t9
      t27 = 0.1q1 / t8
      t28 = 0.1q1 / t9
      t29 = 0.1q1 / t16
      t30 = 0.1q1 / t17
      t11 = 0.1q1 / t11
      t18 = 0.1q1 / t18
      t25 = 0.1q1 / t25
      t26 = 0.1q1 / t26
      t31 = t3 * t25
      t32 = t30 * t26
      t33 = t32 + t31
      t34 = t25 ** 2
      t35 = t25 * t34
      t36 = t3 ** 2
      t37 = t36 ** 2
      t38 = t3 * t36
      t39 = t23 ** 2
      t40 = t1 ** 2
      t41 = t4 ** 2
      t42 = t1 * t12
      t43 = t28 * t27
      t44 = t1 * t11
      t45 = t44 * t23
      t46 = t39 * t11
      t47 = t16 * t11
      t48 = t47 * t23
      t49 = t43 * t1
      t50 = t42 * t18
      t51 = t31 * t23
      t52 = t23 * t26
      t53 = t22 * t28
      t54 = t40 * t27 * t11
      t55 = 0.1q1 / t2
      t20 = 0.1q1 / t20
      t10 = 0.1q1 / t10
      t13 = 0.1q1 / t13
      t56 = 0.1q1 / t14
      t57 = mt ** 2
      t13 = t29 * t55 * t13
      t58 = t22 * t1
      t59 = t40 * t10 * t56
      t60 = t8 * t11
      t61 = t20 ** 2
      t62 = t2 ** 2
      t63 = t26 ** 2
      t64 = t14 ** 2
      t65 = t46 * t43
      t66 = t17 * t55
      t67 = t16 * t63
      t68 = t43 * t11
      t66 = t58 * (t68 * (-t19 * t36 * t17 * t34 + t67) * t61 * t23 * t6
     #4 * t2 + t56 * t21 * t57 * (t10 * (t66 + t60) + (t66 * t27 + t11) 
     #* t28 * t4) + (-t43 * t38 * t17 * t11 * t61 * t34 * t23 + t65 * t3
     #7 * t15 * t17 * t61 * t35) * t64 * t62)
      t69 = t31 * t40
      t70 = t39 * t34
      t63 = t16 ** 2 * t39 * t63
      t71 = t1 * t57
      t19 = t68 * t4 * ((t2 * (t36 * (t24 * t40 * t12 * t34 + t1 * t23 *
     # t25) - t69 * t12 - t70 * t19 * t38) + t62 * (t42 * t23 * t34 * t3
     #8 - t70 * t37) - t40 - t63) * t20 * t14 + t71 * (t23 * (t32 + t31)
     # + t42 * (t13 + t18)))
      t24 = t53 * t20 * t14 * t4 * (t2 * (t45 * t15 * t34 * t36 - t44 * 
     #t31) + t27 * t26 * (t52 * t16 + t1) + t38 * t62 * t23 * t11 * t34)
      t32 = t53 * t27
      t15 = t2 * (t45 * t36 * t17 * t34 * (-t51 * t15 + 1) * t20 * t22 *
     # t14 + t32 * (t26 * (t63 + t40) + t40 * t38 * t15 ** 2 * t17 * t39
     # * t35 + t69 * t17) * t61 * t11 * t64) + t46 * t53 * t64 * t61 * t
     #3 * t37 * t2 * t62 * t17 * t27 * t35 - t46 * t14 * t20 * t37 * t62
     # * t22 * t17 * t35 + t71 * t4 * t55 * t56 * (t43 * t4 + t10)
      t16 = t65 * t21 * t2 * t57 * (t67 * t7 * t30 + t36 * t17 * t34 * (
     #t1 * (t16 * t5 + t3 * t9) - t22 * t3))
      t30 = -4
      t37 = 16
      ret = (t4 * (t2 * (t44 * t28 * (t18 * t27 * t42 + t22 * t25) * t3 
     #- t46 * t34 * t38 + t45 * t25 * (t12 * t25 + t43) * t36) + t49 * (
     #t26 * (-t48 - t22) + t44)) - t53 * t33 * t41 + t54 * (-t1 * t29 + 
     #t52) + t54 * t2 * t28 * (t51 + t50) * t5) * t30 + (t1 * t28 * t56 
     #* (t27 * t3 - t47) * t41 + t59 * t5 + t49 * t11 * (t58 * ((t13 + t
     #18) * t12 * t17 - t29) + (t23 * t33 + t50) * t7 * t2) * t21 * t57 
     #+ t1 * (t56 * (t10 * (-t47 * t8 + t3) + t44) + t43 * (-t48 * t14 *
     # t20 * t26 + t6 * t56)) * t4 + t60 * (-t2 * t38 * t17 * t22 * t39 
     #* t35 + t59) * t9) * t37 - 128 * t66 + 8 * t19 + 24 * t24 - 64 * t
     #15 + 32 * t16 + 48 * t32 * t45 * t57 * t21 * (t31 * t17 + t26) - 1
     #2 * t60 * t36 * t2 * t4 * t22 * t23 * t34

      hjetmass_qqbgg_triangle_pmmp_s124_mhsq_0 = ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     &  hjetmass_qqbgg_triangle_pmmp_s124_mhsq_0_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          real*16 cg

      t1 = za(i2, i3)
      t2 = zb(i2, i1)
      t3 = za(i3, i4)
      t4 = zb(i4, i1)
      t5 = t1 * t2 - t3 * t4
      t6 = za(i1, i3)
      t7 = zb(i3, i1)
      t8 = zb(i3, i2)
      t9 = zb(i4, i3)
      t10 = t6 * t7
      t11 = t1 * t8
      t12 = t3 * t9
      t13 = t10 + t11 + t12
      if ( real(t13) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t11 = cg * sqrt(t13 ** 2) + t10 + t11 + t12
      t12 = za(i1, i2)
      t13 = za(i2, i4)
      t14 = za(i1, i4)
      t15 = zb(i4, i2)
      t16 = t12 * t2
      t17 = t14 * t4
      t18 = t13 * t15 + t16 + t17
      t11 = 0.1q1 / t11
      t19 = -0.2q1 * t1 * t7 * t18 * t11 + t13 * t4
      t20 = 0.2q1 * t3 * t7 * t18 * t11 - t13 * t2
      t10 = -0.2q1 * t10 * t18 * t11 + t16 + t17
      t16 = t6 * (t10 * t7 + t19 * t8)
      t17 = t1 * t15 + t4 * t6
      t21 = 0.2q1 * t6
      t8 = t7 * (t1 * (-t21 * t8 * t18 * t11 + t14 * t15) + t10 * t6)
      t21 = t21 * t18 * t9 * t11 - t12 * t15
      t22 = 0.1q1 / t7
      t18 = 0.1q1 / t18
      t23 = 0.1q1 / t12
      t24 = 0.1q1 / t2
      t14 = 0.1q1 / t14
      t25 = t4 * t23 * t24
      t8 = 0.1q1 / t8
      t21 = 0.1q1 / t21
      t26 = 0.1q1 / t3
      t13 = 0.1q1 / t13
      t16 = 0.1q1 / t16
      t27 = 0.1q1 / t9
      t20 = 0.1q1 / t20
      t28 = t6 * t16
      t21 = t26 * t22 * t21
      t29 = (t20 * t27 + t28) * t19
      t30 = t1 * t10
      t31 = t5 * t1
      t32 = t24 * t23 * t13 * t11
      t33 = t32 * t1
      t34 = t9 * t22
      ret = -0.64q2 * t1 * t4 * t22 * t18 * (t25 + t14) - 0.16q2 * t33 *
     # (t31 * ((-t21 - t8) * t10 * t9 + t26) + (-t30 * t8 - t29) * t17 *
     # t7) - 0.128q3 * t31 * t11 * t18 * (t14 * (t12 * t13 + t34) + (t34
     # * t23 + t13) * t24 * t4) + 0.8q1 * t25 * t13 * t1 * (t30 * (t21 +
     # t8) + t29) + 0.32q2 * t32 * t19 ** 2 * t7 * (t3 * t17 * t20 ** 2 
     #* t27 + t6 ** 2 * t9 * t16 ** 2 * (t1 * (t15 * t3 + t2 * t6) - t5 
     #* t6)) + 0.48q2 * t33 * t5 * t19 * (t28 * t9 + t20)

      hjetmass_qqbgg_triangle_pmmp_s124_mhsq_0_rat = ret/32q0*(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_triangle_pmmp_s12_s123_0
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt
          real*16 cg

      t1 = za(i2, i3)
      t2 = zb(i2, i1)
      t3 = za(i1, i2)
      t4 = za(i3, i4)
      t5 = za(i2, i4)
      t6 = zb(i3, i1)
      t7 = za(i1, i3)
      t8 = zb(i3, i2)
      t9 = t7 * t6
      t10 = t1 * t8
      t11 = t9 + t10
      if ( real(t11) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t11 = cg * sqrt(t11 ** 2) + t10 + t9
      t12 = 0.1q1 / t11
      t13 = 2 * t3
      t14 = t2 * (t13 * t4 * t6 * t12 - t5)
      t15 = za(i1, i4)
      t16 = zb(i4, i1)
      t17 = zb(i4, i2)
      t18 = zb(i4, i3)
      t19 = t3 * t2
      t20 = t19 * (-2 * t9 * t12 + 1)
      t13 = t7 * (-t14 * t18 + t20 * t6) - t10 * t9 * t13 * t2 * t12
      t21 = t1 * t17 + t16 * t7
      t22 = 0.1q1 / t18
      t14 = 0.1q1 / t14
      t11 = 0.1q1 / t11
      t23 = 0.1q1 / t4
      t13 = 0.1q1 / t13
      t24 = t5 * t17
      t25 = t15 * t16
      t26 = t18 ** 2
      t27 = t1 ** 2
      t28 = t1 * t27
      t29 = t13 ** 2
      t30 = t13 * t29
      t31 = t14 ** 2
      t32 = t2 ** 2
      t33 = t7 ** 2
      t34 = t33 ** 2
      t35 = t7 * t33
      t36 = t11 ** 2
      t37 = mt ** 2
      t38 = t12 * (t25 + t24)
      t39 = t2 * t18
      t40 = t33 * t18 * t29
      t41 = t6 * t2
      t42 = t41 * t27
      t43 = t19 + t24 + t25
      t44 = t6 ** 2
      t45 = t19 * t11
      t46 = t43 * t30
      t47 = t10 * t23
      t48 = t12 * t11
      t49 = t48 * t44 * t3
      t50 = 0.1q1 / t3
      t51 = t6 * t17
      t52 = t39 + t51
      t53 = t52 * t12
      t54 = t16 * t18
      t55 = t37 * t16 * t50
      t56 = t23 * t16
      t41 = t27 * (t33 * (-t41 * t16 * t11 * t23 * t13 + t41 * (-t54 * t
     #45 + (t16 * (t53 * t15 - t45 * t8) + t53 * t24 + t53 * t19) * t23 
     #* t1) * t29) + t55 * t23 * t14 * t22 + t23 * (-t1 * t2 * t11 * t52
     # + t55) * t13 * t7 + t56 * t44 * t2 * (t19 * (-t11 + t12) + t38) *
     # t29 * t35)
      t52 = t16 ** 2
      t16 = t8 * t16
      t55 = t3 ** 2
      t57 = t37 * t21 * t12 * t50 * t23
      t53 = t53 * t3
      t5 = t6 * t27 * (t14 * (-t57 * t22 + (-t4 * t6 * t12 * t31 + t11 *
     # t14) * t32 * t1) + t35 * (t48 * t6 * t3 * t32 * ((-t39 - t51 - t1
     #6) * t23 * t1 - t54) * t29 + t23 * t12 * t18 * t6 * t32 * t1 * (-t
     #5 ** 2 * t17 ** 2 - t52 * t15 ** 2 - t32 * t55) * t30) - t56 * t49
     # * t32 * t34 * t29 - t57 * t7 * t13 + t11 * t32 * t1 * (t18 * (t23
     # * t43 - t53) - t53 * t47) * t29 * t33)
      t9 = t12 * t30 * t36 * t18 * t35 * t44 * t55 * t32 ** 2 * t28 * (t
     #9 * t18 + t10 * (t9 * t23 + t18))
      ret = -64 * t42 * (t1 * (t39 * t25 * t24 * t6 * t35 * t12 * t23 * 
     #t30 + t3 * (-t31 * t36 + (t6 * t36 * t23 * t29 + t38 * t6 * t23 * 
     #t30) * t18 * t35 + t33 * t26 * t36 * t29) * t32) + t37 * t6 * t21 
     #* t12 ** 2 * (-t22 * t31 + t40) + t40 * t3 * t32 * t36 * t27 * t8 
     #* t23) + 128 * t49 * t2 * t32 * t28 * (t35 * (-t45 * t4 * t30 * t1
     #8 * t26 + t46 * t26 + t47 * (t19 * (-t10 * t11 + 1) + t24 + t25) *
     # t30 * t18) + t46 * t23 * t6 * t18 * t34 - t45 * t7 * t34 * t44 * 
     #t18 * t23 * t30 - t4 * t14 * t31 * (t45 + 1)) + 16 * t41 - 4 * t1 
     #* (-t14 * (t2 * t27 * t50 * t23 + t22 * t52) + (-t39 * t27 * t50 *
     # t23 - t52 + t56 * (t20 * t50 + t22 * (-t16 + t51)) * t1) * t13 * 
     #t7) + 32 * t5 - 256 * t9 + 8 * t56 * t42 * t13 * t33 * (t13 * t43 
     #- t12)

      hjetmass_qqbgg_triangle_pmmp_s12_s123_0 = ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_triangle_pmmp_s12_s124_0
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt
          real*16 cg

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
      if ( real(t11) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t11 = cg * sqrt(t11 ** 2) + t10 + t9
      t12 = 0.1q1 / t11
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
      t15 = 0.1q1 / t15
      t14 = 0.1q1 / t14
      t23 = 0.1q1 / t18
      t24 = 0.1q1 / t4
      t25 = 0.1q1 / t5
      t26 = t2 ** 2
      t27 = t2 * t26
      t28 = t16 ** 2
      t29 = t7 * t16
      t30 = t8 * t19
      t31 = t13 * t24
      t32 = t1 * t26
      t33 = t2 * t14
      t11 = 0.1q1 / t11
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

      hjetmass_qqbgg_triangle_pmmp_s12_s124_0 = ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_triangle_pmmp_s34_0_0
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt

      t1 = za(i1, i3)
      t2 = za(i2, i4)
      t3 = zb(i3, i1)
      t4 = za(i1, i4)
      t5 = za(i2, i3)
      t6 = zb(i3, i2)
      t7 = t4 * t3
      t8 = t2 * t6
      t9 = zb(i4, i1)
      t10 = zb(i4, i2)
      t11 = za(i1, i2)
      t12 = zb(i2, i1)
      t13 = t5 * t3
      t14 = t6 * t9
      t15 = t3 * t10
      t16 = t3 ** 2
      t17 = -t14 + t15
      t18 = t2 ** 2
      t19 = t8 + t7
      t19 = 0.1q1 / t19
      t20 = 0.1q1 / t11
      t21 = 0.1q1 / t12
      t22 = t19 ** 2
      ret = -16 * za(i3, i4) * (-t1 ** 2 * t2 * t3 * t16 + t11 * t12 * (
     #t2 * (-2 * t1 * t16 + t2 * (-t14 + 2 * t15) - t13 * t6) + t7 * (t2
     # * t9 + t13)) + t10 * t17 * t2 * t18 + t1 * t5 * t16 * (-t8 + t7) 
     #+ t4 * t9 * t17 * t18 + t4 * t5 ** 2 * t16 * t6) * zb(i4, i3) * t2
     #0 * t21 * t19 * t22

      hjetmass_qqbgg_triangle_pmmp_s34_0_0 = ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_triangle_pmmp_s34_mhsq_s12
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt

      t1 = za(i2, i3)
      t2 = za(i1, i3)
      t3 = zb(i3, i1)
      t4 = zb(i3, i2)
      t5 = za(i1, i4)
      t6 = zb(i4, i1)
      t7 = za(i2, i4)
      t8 = zb(i4, i2)
      t9 = t2 * t3
      t10 = t1 * t4
      t11 = t5 * t6
      t12 = t7 * t8
      t13 = t11 + t12 + t9 + t10
      t14 = za(i1, i2)
      t15 = za(i3, i4)
      t16 = zb(i2, i1)
      t17 = zb(i4, i3)
      t13 = -4 * t14 * t15 * t16 * t17 + t13 ** 2
      t13 = t11 + t12 + sqrt(t13) + t9 + t10
      t18 = (0.1q1 / 0.4q1)
      t19 = t13 ** 2
      t20 = t14 * t15
      t21 = t20 * t16 * t17
      t22 = t19 * t18 - t21
      t22 = 0.1q1 / t22
      t23 = (0.1q1 / 0.2q1)
      t24 = t13 * t18
      t25 = t23 * t2 * t16 * t17 - t24 * t8
      t26 = t13 * t14 * t22
      t27 = t26 * t25
      t28 = t23 * t1 * t16 * t17 + t24 * t6
      t29 = t13 * t15 * t22
      t30 = t29 * t28
      t31 = t1 * t3 + t6 * t7
      t19 = t19 * t22
      t32 = t19 * t31
      t33 = t2 * t4 + t5 * t8
      t34 = t26 * t16
      t35 = t11 + t9
      t36 = t19 * t18
      t37 = t36 * t14 * t16
      t38 = -t23 * t34 * t35 + t37
      t39 = t21 * t23 * t13 * t22
      t35 = t18 * t19 * t35 - t39
      t40 = t38 * t35
      t41 = -t32 * t34 * t33 / 8
      t42 = -t30 * t26 * (t23 * t5 * t16 * t17 + t24 * t4) + t40 + t41
      t41 = t27 * t29 * (t23 * t7 * t16 * t17 - t24 * t3) + t40 + t41
      t43 = 2 * t14 * t16 + t13
      t21 = t23 * t13 ** 2 - 2 * t21
      t44 = t12 + t10
      t37 = -t23 * t34 * t44 + t37
      t9 = t9 + t10
      t10 = t18 * t19 * t9 - t39
      t45 = 2 * t15 * t17 + t13
      t31 = t34 * t31
      t20 = t20 * t23
      t34 = t24 * t1 + t20 * t6
      t46 = t13 * t16 * t22
      t47 = -t46 * t34
      t48 = -t24 * t2 + t20 * t8
      t22 = t13 * t17 * t22
      t49 = t22 * t48
      t2 = t1 * t8 + t2 * t6
      t8 = t29 * t17
      t50 = t8 * t2
      t33 = -t31 * t19 * t33 / 8
      t4 = t47 * t22 * (t20 * t4 + t24 * t5) + t40 + t33
      t5 = t26 * t28
      t3 = t49 * t46 * (t20 * t3 - t24 * t7) + t40 + t33
      t7 = t18 * t19 * t44 - t39
      t11 = t11 + t12
      t12 = t36 * t15 * t17
      t20 = -t23 * t8 * t11 + t12
      t2 = t19 * t2
      t22 = -t22 * t34
      t8 = -t23 * t8 * t9 + t12
      t9 = t18 * t19 * t11 - t39
      t11 = -t29 * t25
      t4 = 0.1q1 / t4
      t12 = 0.1q1 / t41
      t18 = 0.1q1 / t21
      t3 = 0.1q1 / t3
      t19 = 0.1q1 / t17
      t21 = 0.1q1 / t42
      t24 = 0.1q1 / t15
      t25 = 0.1q1 / t14
      t26 = 0.1q1 / t16
      t28 = t21 ** 2
      t29 = t21 * t28
      t33 = t12 ** 2
      t34 = t12 * t33
      t36 = t28 - t33
      t39 = t7 ** 2
      t41 = t20 ** 2 + t39
      t42 = t37 + t9
      t44 = -t7 - t8
      t51 = t4 - t3
      t52 = t35 + t7
      t53 = t38 + t37
      t54 = t35 ** 2
      t55 = t35 * t54
      t56 = t27 ** 2
      t57 = t30 ** 2
      t58 = t3 ** 2
      t59 = t3 * t58
      t60 = t1 ** 2
      t61 = t8 ** 2
      t62 = t6 ** 2
      t63 = t4 ** 2
      t64 = t4 * t63
      t65 = t18 ** 2
      t66 = t45 ** 2
      t67 = t20 * t2
      t68 = t50 * t53
      t69 = t49 * t47
      t70 = t30 * t27
      t71 = t70 * t32
      t72 = t71 * t38
      t73 = t72 * t5 * t11
      t74 = t73 * t33
      t75 = t69 * t35
      t76 = t75 * t31
      t77 = t76 * t50 * t58
      t78 = t56 * t57 * t32
      t79 = t24 * t19
      t80 = t79 * t16 * t14
      t81 = t52 * t33
      t82 = t51 * t25 * t24
      t13 = 0.1q1 / t13
      t83 = t35 + t7 + t8
      t84 = t38 + t37 + t10
      t85 = -t35 - t7 - t20
      t86 = t38 ** 2
      t87 = t10 ** 2
      t88 = mt ** 2
      t89 = t30 * t38
      t90 = t35 * t28
      t91 = t38 * (t29 - t34)
      t92 = t85 * t33
      t93 = t83 * t28
      t94 = t18 * t45
      t95 = t18 * t47
      t96 = t95 * t43
      t97 = t80 * t30
      t15 = t15 * t17
      t17 = t15 * t43 ** 2
      t98 = t17 * t65
      t14 = t14 * t16
      t16 = t79 * t35
      t99 = t89 * t21
      t100 = t35 * t47
      t101 = t100 * t4
      t102 = t38 + t37 + t9
      t103 = t37 + t10
      t104 = t9 ** 2
      t105 = t20 * t58
      t106 = t44 * t63
      t107 = t94 * t47
      t108 = t107 * t54 * t63
      t109 = t30 * t4
      t110 = t109 * t94
      t111 = t30 * t2
      t112 = t50 * t9
      t113 = t5 * t28
      t114 = t42 * t33
      t115 = t28 * t103
      t116 = t5 * t12
      t117 = t98 * t26 * t25
      t118 = t117 * t30 * (-t116 * t38 * t2 + (t102 * t33 * t47 + t89 * 
     #(-t104 * t34 + t29 * t86)) * t32 * t56)
      t11 = t45 * (t38 * (t111 * t22 * t12 - t77) + t78 * (-t114 + t115)
     # + t72 * (t113 * t11 + t112 * t33)) * t65 * t43
      t119 = t94 * t80 * (-t94 * t78 * t38 * t54 * t34 + (t100 * (t94 * 
     #(t105 + t106) + t63) - t108 + t110) * t50 * t49 * t31)
      t120 = t70 * t60 * t25 * t24
      t121 = t37 ** 2
      t112 = (-t69 * t62 * t3 + t70 * (t12 * (t112 * t88 * t32 * t38 * t
     #13 * t24 * t25 * t12 - t62) + t21 * t62)) * t26 * t19
      t122 = t117 * (-t100 * t5 * t2 * t3 + t91 * t78 * t121)
      t11 = t45 * (t80 * (t78 * t36 - t77) * t18 + t43 * (t69 * (t58 * (
     #-t42 * t50 + t67) + t68 * t63) * t31 * t35 - t74) * t65) + t69 * (
     #t62 * t26 * t19 * t4 + t82 * t60) + t80 * (t69 * t31 * t54 * t50 *
     # t58 + t78 * (-t38 * t41 * t34 + t28 * t44 + t81 + (t54 + t39 + t6
     #1) * t38 * t29)) * t65 * t66 + t98 * t2 * t5 * t25 * t26 * (t99 + 
     #t101) + (-t97 * t66 * t50 * t65 * t3 + t16 * t47 * (t67 * t88 * t1
     #3 * t25 * t26 + t14 * t66 * t50 * t7 * t65) * t58) * t49 * t31 + t
     #30 * (t96 * (-t28 + t33 + t94 * (t92 + t93)) + t97 * (t91 + (t20 *
     # t33 - t90) * t65 * t66) + t98 * (-t84 * t28 * t47 + t89 * (t29 * 
     #t87 - t34 * t86)) * t26 * t25) * t32 * t56 + t120 * t21 + t11 + t1
     #18 + t119 - t99 * t43 * t65 * t45 * t2 * t22 + t122 + t112 - t120 
     #* t12
      t44 = t49 ** 2
      t56 = t15 * t43
      t60 = t56 * t37 * t25 * t26
      t62 = t94 * t35
      t67 = t80 * t45
      t77 = t18 * t43
      t91 = t43 * t63
      t46 = -t22 * t46 * t48
      t48 = t45 * t65
      t97 = t78 * t18
      t112 = t97 * t38
      t118 = -t59 + t64
      t119 = -t7 - t20
      t120 = t28 * t38
      t122 = t78 * t38
      t123 = t79 * t88
      t124 = t123 * t13
      t125 = t124 * t25 * t26
      t126 = t79 * t18
      t127 = t48 * t43
      t39 = -t69 * t50 * t3 * (t125 + t127) + t95 * (t14 * (t79 * t30 * 
     #t63 * t45 + t126 * (-t30 * t8 * t63 - t55 * t47 * t59 + t100 * (t3
     #9 * t64 - t41 * t59)) * t66) + t17 * t95 * t26 * t25 * (t38 * t58 
     #+ (t87 + t86 + t121) * t64 * t35)) * t44
      t41 = t88 * t13
      t87 = t125 * t2
      t24 = t78 * (t38 * (t19 * (-t41 * t24 * t25 * t26 * t28 + t94 * t1
     #4 * (t35 + t8) * t24 * t29) + t77 * t33 * (t12 * t42 + t94)) + t77
     # * (t34 * (t56 * t9 * t18 * t25 * t26 + 1) - t29) * t86)
      t128 = t100 * (-t48 * t51 * t5 * t50 * t43 - t87 * t22 * t3)
      t24 = t31 * t39 + t31 * (t48 * t69 * t43 * t50 * t4 + t47 * (t80 *
     # (-t63 * t52 * t65 * t66 * t30 + t100 * t118) + t65 * t47 * t43 * 
     #(t45 * (t119 * t58 - t106) - t91 * t15 * t84 * t26 * t25)) * t44) 
     #+ t43 * (t48 * (t100 * t2 * t22 * t51 + t78 * (t86 * (t119 * t34 +
     # t29 * t83) - t120)) - t112 * t29 * t103) - t79 * t73 * t88 * t13 
     #* t25 * t26 * t28 + t48 * (t31 * (t47 * t30 * (-t67 * t85 * t58 + 
     #t91 * t53) * t44 + t91 * t75 * (-t10 * t50 + t46)) + t67 * t22 * t
     #50 * (-t100 * t3 + t99)) + t112 * (t67 * (t29 * t7 + t34 * (t7 * (
     #-1 + t62) - t20 - t35)) + t77 * (t34 * (t9 * (-t20 * t45 + t60) - 
     #t40 * t45) + t45 * (t10 * t83 + t37 * t83) * t29)) + t125 * (t101 
     #* t2 * t22 + t122 * t33) + t24 + t128
      t39 = t80 * t66
      t40 = t39 * t65
      t73 = -t37 - t10
      t75 = t47 ** 2
      t78 = t39 * t35
      t106 = t12 - t21
      t112 = t17 * t5 * t25 * t26
      t128 = t77 * t31
      t129 = t38 * t22
      t130 = t31 * t54
      t131 = t38 * t57
      t132 = t131 * t27
      t133 = t63 - t58
      t134 = t35 * t22
      t135 = t132 * t22
      t136 = t69 * t30
      t41 = t41 * t35
      t137 = t35 * t5
      t81 = t79 * (-t41 * t31 * t75 * t44 * t58 * t25 * t26 + t14 * (t94
     # * (-t130 * t118 * t44 * t75 - t134 * t133 * t49 * t75 - t135 * t2
     #8) + (t75 * (t130 * (t119 * t59 + t64 * t7) * t44 - t134 * t58 * t
     #52 * t49) + t136 * t22 * t3 + t135 * (t90 - t81)) * t65 * t66))
      t62 = t77 * (t27 * (t95 * (t56 * t106 * t26 * t25 * t5 - t45 * t22
     # * t12) * t30 + t131 * (t33 * (t129 * t94 - t5) + t113)) + t62 * t
     #49 * t75 * (t133 * t49 * t31 + t137 * t58))
      t90 = -t38 - t37 - t9
      t118 = t7 + t8
      t119 = t4 * t35
      t131 = t31 * t35
      t134 = (t37 + t9) * t33
      t136 = t75 * (t43 * (t56 * t18 * t4 * t5 * (t119 * (t38 + t37 + t1
     #0) - 1) * t26 * t25 + t94 * (t4 - t3) * t22 + t137 * (t63 - t58)) 
     #* t49 + t131 * (t43 * (t38 * t64 + t59 * t90) + (t80 * t59 - t77 *
     # t64 * (t38 + t37)) * t7 * t45) * t44) - t136 * t94 * t43 * t5 * t
     #3 + t18 * t30 * (-t45 * t43 * t38 * t5 * t50 * t21 + t70 * (t43 * 
     #(t45 * (t21 * (-t86 * t22 * t21 + t5) - t116) - t134 * t56 * t38 *
     # t5 * t26 * t25) + t129 * t39 * t28 * t118))
      t137 = t69 * t4
      t138 = t38 * t50
      t139 = t49 * t31
      t140 = t49 * t75
      t141 = t20 * t59
      t142 = t8 * t64
      t143 = t38 * t37
      t105 = t94 * (t77 * (t30 * (t5 * (t138 * t12 + t137) + t27 * t47 *
     # t22 * t21) + t140 * (t139 * t9 * t59 - t5 * t63) * t54 + t140 * (
     #t5 * (-t63 * t7 + t105) + t58 * t42 * t22) * t35) + t80 * (t75 * (
     #t94 * t54 * t22 * t63 * t49 + t131 * (t141 + t94 * (-t141 + t142) 
     #* t7) * t44) - t110 * t69 * t22 + t135 * t33))
      t36 = (-t124 * t99 * t50 * t5 + t98 * (t57 * t27 * t86 * t5 * t36 
     #+ t131 * (t59 * (-t53 * t9 - t143) + t64 * (t10 * t53 + t143)) * t
     #44 * t75)) * t26 * t25
      t14 = t140 * t18 * t35 * (t112 * t58 * t18 * t90 + t14 * (t66 * (t
     #126 * t22 * t118 * t63 + t142 * t139 * t16 * t18) - t139 * t79 * t
     #118 * t64 * t45))
      t1 = t14 + t18 * t136 + t75 * (t128 * t35 * (t103 * t64 + t94 * (t
     #59 * (t102 * t20 + t102 * t7) + t64 * (-t10 * t7 - t8 * t84))) * t
     #44 + t65 * t3 * (-t78 * t20 * t22 * t3 + t112) * t49) + t123 * t69
     # * t26 * t25 * (t46 * t35 * t31 * t63 * t13 - t51 * t6 * t1) + t70
     # * (t123 * t1 * t6 * t25 * t26 * t106 + (t113 * t17 * t38 * t25 * 
     #t26 * t103 + t39 * (t12 * (-t38 * t20 * t12 + 1) - t21) * t22) * t
     #65 * t30) + t125 * (t76 * (-t46 * t58 + t69 * t63) + t116 * t89 * 
     #t50) + t127 * (t75 * (t35 * (t58 * (t5 * t7 + t129) + t63 * (-t22 
     #* t84 - t5 * t8)) * t49 + t130 * (t53 * t59 - t64 * t84) * t44) + 
     #t132 * (t22 * (t114 - t115) + t5 * (-t92 - t93))) + t81 + t62 + t3
     #6 + t105
      t5 = t40 * t50
      t6 = t111 * t106
      t5 = t71 * (t45 * (t138 * t80 * t28 * t18 + t43 * t38 * (-t8 * t2 
     #* t28 - t68 * t33) * t65) + t5 * (t21 * (-t38 * t21 * t83 + 1) - t
     #12) - t120 * t87 * t8) + t2 * t31 * (t41 * t82 * t50 * t26 * t19 +
     # t137 * t117 * (-t119 * t84 + 1) + (-t100 * t63 * t18 + t48 * (t10
     #0 * t7 * t133 - t47 * t54 * t58 - t109)) * t49 * t43) + t128 * t2 
     #* t49 * (t15 * t96 * t3 * (t3 * t35 * t102 - 1) * t26 * t25 + t108
     # + t94 * t30 * t3 + t100 * (t94 * t8 * t63 + t58)) + (t45 * (-t89 
     #* t80 * t33 * t50 * t18 + t43 * (t50 * (t106 * t47 + t30 * (t120 *
     # t103 + t28 * t86)) - t6) * t65) - t6 * t125 - t92 * t5 * t89) * t
     #32 * t27
      t6 = t33 - t28
      t14 = -t12 + t21
      t3 = t2 * (t43 * (t72 * t6 * t18 + t48 * (t72 * (t28 * (t35 + t7) 
     #+ t33 * (-t35 - t7 - t20)) + t131 * t50 * (-t4 + t3))) + t117 * t3
     #2 * t27 * (t14 * t47 + t30 * (t38 * (t28 * t73 + t134) + t6 * t86)
     #) + t138 * t125 * t32 * t14)
      t6 = 8
      ret = t23 * t138 * t127 * t32 * t2 * t106 + t6 * (t18 * (t31 * (t9
     #1 * t107 * t30 * t10 * t44 + (t43 * t58 + t63 * (t78 * t61 * t4 * 
     #t18 - t43) + t17 * t59 * t18 * t35 * (-t104 - t86 - t121) * t26 * 
     #t25) * t44 * t75) + t97 * t17 * t26 * t25 * t86 * (t29 * t73 + t34
     # * t37)) + t31 * (t95 * (t67 * (t107 * t55 * t64 - t30 * t58) + t7
     #7 * t58 * (t56 * t42 * t26 * t25 * t47 - t45 * t30 * t102)) * t44 
     #+ t69 * (-t46 * t127 * t35 * t58 + (-t16 * t10 * t13 * t25 * t26 *
     # t63 + t79 * t26 * t25 * t13 * t4) * t50 * t88)) + t24 + t125 * t7
     #4 + t122 * t65 * (t43 * (-t60 * t10 * t29 + (t37 * t85 - t52 * t9)
     # * t34 * t45) + t39 * (t29 * (-t35 * t7 - t52 * t8) + t52 * t34 * 
     #t20)) + t40 * t22 * t50 * (-t89 * t12 + t101)) + 16 * t1 - 4 * t11
     # + 2 * t5 - t3

      hjetmass_qqbgg_triangle_pmmp_s34_mhsq_s12 = ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     &  hjetmass_qqbgg_triangle_pmmp_s34_mhsq_s12_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          real*16 cg

          cg = 1q0

      t1 = za(i1, i3)
      t2 = zb(i3, i1)
      t3 = za(i2, i3)
      t4 = zb(i3, i2)
      t5 = za(i1, i4)
      t6 = zb(i4, i1)
      t7 = za(i2, i4)
      t8 = zb(i4, i2)
      t9 = t1 * t2
      t10 = t3 * t4
      t11 = t5 * t6
      t12 = t7 * t8
      t13 = t11 + t12 + t9 + t10
      t14 = za(i1, i2)
      t15 = za(i3, i4)
      t16 = zb(i2, i1)
      t17 = zb(i4, i3)
      t18 = t14 * t16
      t13 = -0.4q1 * t18 * t15 * t17 + t13 ** 2
      t13 = cg * sqrt(t13) + t10 + t11 + t12 + t9
      t19 = t3 * t2
      t20 = t7 * t6
      t21 = t20 + t19
      t22 = t18 * t21
      t23 = t1 * t6 + t3 * t8
      t24 = t11 + t9
      t25 = 0.1q1 / 0.2q1
      t26 = t18 * t15 * t17
      t27 = t13 * t24 * t25 - t26
      t28 = 0.1q1 / 0.4q1
      t29 = t13 ** 2
      t30 = t28 * t29 - t26
      t31 = t11 + t9
      t32 = t15 ** 2
      t33 = t31 * t16
      t34 = t33 * t17 * t15
      t35 = t4 * t6
      t36 = t11 * t9 * t13
      t19 = t20 + t19
      t20 = t1 ** 2 * t2 ** 2
      t37 = t5 ** 2 * t6 ** 2
      t4 = t4 * t19 * t1
      t19 = t8 * t19 * t5
      t38 = t18 * t13
      t31 = t31 * t14
      t5 = t3 * t5
      t39 = t29 * t16
      t40 = t25 * t38 * ((t11 + t18 + t10) * t17 * t15 + t20 + t37 + t4 
     #+ t19)
      t41 = -t28 * t39 * (-t5 * t17 + t31) - t18 * (t14 * (-t35 * t17 * 
     #t32 + t34) - t36) + t40
      t2 = t2 * t8
      t7 = t1 * t7
      t4 = t25 * t38 * ((t12 + t18 + t9) * t17 * t15 + t20 + t37 + t4 + 
     #t19)
      t19 = -t28 * t39 * (t7 * t17 + t31) - t18 * (t14 * (t2 * t17 * t32
     # + t34) - t36) + t4
      t20 = -t18 * t24 + t38 * t25
      t24 = t17 ** 2
      t31 = t31 * t17
      t32 = t29 * t14
      t5 = t28 * t32 * (t16 * (-t11 - t9) + t35 * t15) + t40 + t18 * ((t
     #5 * t24 - t31) * t16 * t15 + t36)
      t2 = -t28 * t32 * (t2 * t15 + t33) + t4 - t18 * ((t7 * t24 + t31) 
     #* t16 * t15 - t36)
      t4 = t14 * t15
      t7 = t25 * t13
      t18 = t7 * t3 + t4 * t6
      t24 = -t16 * t18
      t4 = -t1 * t7 + t4 * t8
      t28 = t17 * t4
      t1 = t1 * t16 * t17 - t7 * t8
      t8 = t14 * t1
      t31 = t3 * t16 * t17 + t7 * t6
      t32 = t15 * t31
      t18 = -t17 * t18
      t33 = t15 * t17
      t34 = t33 * t23
      t31 = t31 * t14
      t11 = t11 + t12
      t7 = t33 * t7
      t9 = t9 + t10
      t10 = 0.1q1 / t17
      t5 = 0.1q1 / t5
      t12 = 0.1q1 / t14
      t14 = 0.1q1 / t15
      t17 = 0.1q1 / t19
      t2 = 0.1q1 / t2
      t19 = 0.1q1 / t41
      t35 = 0.1q1 / t16
      t36 = t2 ** 2
      t37 = t5 ** 2 - t36
      t30 = 0.1q1 / t30 ** 2
      t12 = t35 * t12
      t30 = t12 * t30
      t35 = t8 * t32
      t38 = t5 - t2
      t39 = t27 * (-t17 ** 2 + t19 ** 2)
      ret = -0.16q2 * t12 * t14 * t10 * t6 * t3 * ((t17 - t19) * t28 * t
     #24 + t35 * (-t5 + t2)) + 0.8q1 * t30 * t14 * t10 * t13 * (t22 * (t
     #39 * t28 ** 2 * t24 ** 2 + (t34 * (t19 * (-t27 * (t25 * t13 * t9 -
     # t26) * t19 - 0.1q1) + t17) - t39 * t16 * t4 * t18) * t28 * t24) +
     # t20 * t31 * t32 * t34 * t38) + 0.2q1 * t30 * t21 * t23 * t13 * t2
     #9 * (t20 * t23 * t38 + t35 * (t5 * (-t20 * (-t33 * t9 + t7) * t5 -
     # 0.1q1) + t2) * t14 * t10) - 0.4q1 * t30 * t29 * (((t37 * t32 ** 2
     # * t8 ** 2 + (-t37 * t15 * t1 * t31 + t34 * (t25 * t13 * t11 - t26
     #) * t36) * t32 * t8) * t20 * t21 + t24 * t27 * t23 * (t17 * (t22 *
     # t28 * (-t33 * t11 + t7) * t17 - t18) + t18 * t19)) * t14 * t10 + 
     #t22 * t23 ** 2 * t27 * (-t17 + t19))

      hjetmass_qqbgg_triangle_pmmp_s34_mhsq_s12_rat = ret/32q0*(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_triangle_pmpm_s123_mhsq_0
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt
          real*16 cg

      t1 = za(i1, i4)
      t2 = za(i2, i4)
      t3 = zb(i3, i1)
      t4 = zb(i4, i1)
      t5 = za(i1, i2)
      t6 = zb(i2, i1)
      t7 = za(i1, i3)
      t8 = za(i2, i3)
      t9 = zb(i3, i2)
      t10 = t5 * t6
      t11 = t7 * t3
      t12 = t8 * t9 + t10 + t11
      t13 = zb(i4, i2)
      t14 = za(i3, i4)
      t15 = zb(i4, i3)
      t16 = t1 * t4
      t17 = t2 * t13
      t18 = t14 * t15
      t19 = t16 + t17 + t18
      if ( real(t19) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t18 = cg * sqrt(t19 ** 2) + t16 + t17 + t18
      t19 = 0.1q1 / t18
      t11 = -2 * t16 * t12 * t19 + t10 + t11
      t20 = -2 * t1 * t12
      t7 = t4 * (t1 * t11 + t2 * (t20 * t13 * t19 + t7 * t9))
      t21 = -2 * t2 * t12 * t4 * t19 + t3 * t8
      t22 = -2 * t14 * t12 * t4 * t19 - t6 * t8
      t23 = t14 * t3 + t2 * t6
      t24 = t21 * t13
      t25 = t1 * (t11 * t4 + t24)
      t20 = t20 * t15 * t19 - t5 * t9
      t26 = t1 * t3 + t2 * t9
      t8 = 0.1q1 / t8
      t18 = 0.1q1 / t18
      t27 = 0.1q1 / t6
      t28 = 0.1q1 / t12
      t29 = 0.1q1 / t4
      t22 = 0.1q1 / t22
      t25 = 0.1q1 / t25
      t30 = 0.1q1 / t5
      t31 = t3 ** 2
      t32 = t18 ** 2
      t33 = t25 ** 2
      t34 = t25 * t33
      t35 = t2 ** 2
      t36 = t22 ** 2
      t37 = t21 ** 2
      t38 = t1 ** 2
      t39 = t38 ** 2
      t40 = t1 * t38
      t41 = t12 ** 2
      t42 = t4 ** 2
      t43 = mt ** 2
      t44 = t14 ** 2 * t37 * t36
      t45 = t1 * t25
      t46 = t30 * t27
      t47 = t46 * t41 * t32
      t48 = t2 * t8
      t49 = t48 * t38
      t50 = t46 * t43
      t51 = t18 * t12
      t52 = t40 * t15 * t8
      t36 = t14 * t36
      t53 = t8 * t30
      t54 = t27 * t23
      t55 = 0.1q1 / t14
      t7 = 0.1q1 / t7
      t56 = 0.1q1 / t15
      t20 = 0.1q1 / t20
      t57 = t2 * t21
      t58 = t11 * t25
      t59 = t56 * t22
      t20 = t55 * t29 * t20
      t60 = t2 * t11
      t61 = t46 * t8
      t62 = t61 * t3
      t57 = t62 * (t51 * (t40 * (-t57 * t11 * t33 * t42 + t17 * t37 * t3
     #3 * t4) + t35 + t44 + t39 * t42 * t37 * t33 + (-t24 * t35 * t11 * 
     #t33 - t57 * t25) * t4 * t38 + t58 * t16 * t35) + t43 * t2 * (t21 *
     # (t59 - t45) + t60 * (t20 - t7)))
      t59 = t59 - t45
      t60 = t60 * t7
      t63 = t14 * t8
      t64 = t21 * t22
      t65 = t50 * t48
      t7 = t2 * (-t27 * t28 * (t1 * t30 + t63) * t31 - t64 * t62 * t51 *
     # t14) - t28 * (t46 * t9 + t8) * t3 * t35 + t65 * ((t21 * t59 - t60
     #) * t26 * t4 + t23 * t2 * (-t55 + (-t20 + t7) * t11 * t15)) * t19 
     #- t10 * t52 * t23 * t4 * t37 * t34
      t10 = t48 * t16
      t11 = t54 * t51 * t3 * (t24 * t49 * t4 * t33 + t30 * t22 * (-t64 *
     # t14 + t2) + t40 * t42 * t21 * t8 * t33 - t10 * t25)
      t10 = t3 * (t46 * t2 * (t22 * (t63 * t21 - t23) + t48) - t40 * t4 
     #* t37 * t8 * t33 + t10 * t27 * (t23 * t25 + t60 * t30) + t49 * t21
     # * t25 * (t58 + t46) * t4) + t53 * t35 * (t2 * t55 + t64) + t61 * 
     #t35 * t4 * (t45 * t21 + t60) * t9 + t54 * t59 * t31
      t16 = 64
      ret = t16 * (t4 * (t49 * t15 * t21 * t33 * (t24 * t45 - 1) * t18 *
     # t23 * t12 + t47 * (t22 * (t35 + t44) - t40 * t35 * t13 ** 2 * t15
     # * t37 * t34 - t45 * t35 * t15) * t8 * t23) - t47 * t1 * t39 * t23
     # * t4 * t42 * t37 * t15 * t8 * t34 + t50 * t2 * t31 * t28 * t29 + 
     #t51 * t39 * t23 * t42 * t37 * t15 * t8 * t34) - 128 * t54 * t2 * (
     #t53 * (-t17 * t38 * t15 * t33 + t36) * t32 * t21 * t41 * t4 + (t39
     # * t13 * t15 * t8 * t32 * t34 * t30 * t37 - t52 * t32 * t33 * t30 
     #* t21) * t41 * t42 + t43 * t3 * t19 * t28 * (t15 * t30 * t29 - t8)
     #) + 8 * t57 + 16 * t7 - 32 * t50 * t19 * t8 * t37 * t4 * (t36 * t2
     #6 * t56 + t38 * t15 * t33 * (t1 * t23 - t2 * (t1 * t6 - t14 * t9))
     #) + 48 * t65 * t23 * t21 * t19 * (t45 * t15 - t22) - 24 * t11 + 4 
     #* t10 + 12 * t5 * t38 * t3 * t23 * t4 * t21 * t8 * t33

      hjetmass_qqbgg_triangle_pmpm_s123_mhsq_0 = ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     &  hjetmass_qqbgg_triangle_pmpm_s123_mhsq_0_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          real*16 cg

      t1 = za(i2, i4)
      t2 = zb(i2, i1)
      t3 = za(i3, i4)
      t4 = zb(i3, i1)
      t5 = t1 * t2 + t3 * t4
      t6 = za(i1, i4)
      t7 = zb(i4, i1)
      t8 = zb(i4, i2)
      t9 = zb(i4, i3)
      t10 = t6 * t7
      t11 = t1 * t8
      t12 = t3 * t9
      t13 = t10 + t11 + t12
      if ( real(t13) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t11 = cg * sqrt(t13 ** 2) + t10 + t11 + t12
      t12 = za(i1, i2)
      t13 = za(i2, i3)
      t14 = za(i1, i3)
      t15 = zb(i3, i2)
      t16 = t12 * t2
      t17 = t14 * t4
      t18 = t13 * t15 + t16 + t17
      t11 = 0.1q1 / t11
      t19 = -0.2q1 * t1 * t18 * t7 * t11 + t13 * t4
      t20 = -0.2q1 * t3 * t18 * t7 * t11 - t13 * t2
      t10 = -0.2q1 * t10 * t18 * t11 + t16 + t17
      t16 = t6 * (t10 * t7 + t19 * t8)
      t17 = t1 * t15 + t4 * t6
      t21 = -0.2q1 * t6 * t18
      t8 = t7 * (t1 * (t21 * t8 * t11 + t14 * t15) + t10 * t6)
      t14 = t21 * t9 * t11 - t12 * t15
      t20 = 0.1q1 / t20
      t12 = 0.1q1 / t12
      t21 = 0.1q1 / t2
      t16 = 0.1q1 / t16
      t13 = 0.1q1 / t13
      t22 = t6 * t16
      t23 = t5 * t1
      t24 = 0.1q1 / t3
      t8 = 0.1q1 / t8
      t14 = 0.1q1 / t14
      t25 = 0.1q1 / t7
      t26 = 0.1q1 / t9
      t14 = t24 * t25 * t14
      t27 = (-t20 * t26 + t22) * t19
      t28 = t1 * t10
      t29 = t13 * t21 * t12
      t30 = t29 * t11
      t18 = 0.1q1 / t18
      ret = -0.48q2 * t23 * t19 * t11 * t12 * t13 * t21 * (-t22 * t9 + t
     #20) - 0.8q1 * t29 * t4 * t1 * (t28 * (-t14 + t8) + t27) - 0.16q2 *
     # t30 * t1 * (t23 * ((t14 - t8) * t10 * t9 + t24) + (t28 * t8 + t27
     #) * t17 * t7) + 0.128q3 * t23 * t4 * t11 * t21 * t18 * (-t9 * t12 
     #* t25 + t13) + 0.32q2 * t30 * t19 ** 2 * t7 * (-t3 * t17 * t20 ** 
     #2 * t26 + t6 ** 2 * t9 * t16 ** 2 * (t1 * (-t15 * t3 + t2 * t6) - 
     #t5 * t6)) + 0.64q2 * t1 * t4 ** 2 * t12 * t21 * t18 * t25

      hjetmass_qqbgg_triangle_pmpm_s123_mhsq_0_rat = ret/32q0*(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_triangle_pmpm_s124_mhsq_0
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt
          real*16 cg

      t1 = za(i2, i4)
      t2 = zb(i3, i1)
      t3 = za(i1, i2)
      t4 = zb(i3, i2)
      t5 = za(i1, i4)
      t6 = zb(i4, i3)
      t7 = t3 * t4 - t5 * t6
      t8 = za(i1, i3)
      t9 = zb(i2, i1)
      t10 = zb(i4, i1)
      t11 = zb(i4, i2)
      t12 = t3 * t9
      t13 = t5 * t10
      t14 = t1 * t11 + t12 + t13
      t15 = za(i2, i3)
      t16 = za(i3, i4)
      t17 = t8 * t2
      t18 = t15 * t4
      t19 = t16 * t6
      t20 = t17 + t18 + t19
      if ( real(t20) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t19 = cg * sqrt(t20 ** 2) + t17 + t18 + t19
      t20 = 0.1q1 / t19
      t13 = -2 * t17 * t14 * t20 + t12 + t13
      t21 = 2 * t8
      t22 = -t21 * t4 * t14 * t20 + t11 * t5
      t23 = t8 * t13
      t24 = t2 * (t15 * t22 + t23)
      t11 = t21 * t14 * t6 * t20 - t11 * t3
      t21 = t1 * t6
      t25 = t3 * t2 + t21
      t26 = -2 * t15 * t2 * t14 * t20 + t1 * t10
      t27 = 2 * t16 * t2 * t14 * t20 - t1 * t9
      t28 = t2 * t13
      t29 = t8 * (t26 * t4 + t28)
      t5 = t5 * t2
      t30 = t1 * t4
      t31 = t30 + t5
      t32 = 0.1q1 / t6
      t10 = 0.1q1 / t10
      t33 = 0.1q1 / t9
      t11 = 0.1q1 / t11
      t3 = 0.1q1 / t3
      t34 = 0.1q1 / t15
      t19 = 0.1q1 / t19
      t27 = 0.1q1 / t27
      t35 = 0.1q1 / t16
      t36 = 0.1q1 / t14
      t24 = 0.1q1 / t24
      t29 = 0.1q1 / t29
      t37 = t2 * t24
      t38 = t2 ** 2
      t39 = t2 * t38
      t40 = t13 ** 2
      t41 = t24 ** 2
      t42 = t24 * t41
      t43 = mt ** 2
      t44 = t13 * (t11 * t35 + t37)
      t45 = t2 * t32
      t46 = t21 * t3
      t47 = t3 * t33
      t48 = t47 * t10
      t49 = t1 * t38
      t50 = t33 * t11
      t51 = t3 * t1
      t52 = t8 * t29
      t53 = t38 * t33 * t10
      t54 = t19 ** 2
      t55 = t15 ** 2
      t56 = t11 ** 2
      t57 = t14 ** 2
      t58 = t6 ** 2 * t40 * t56
      t59 = t39 * t16
      t60 = t59 * t4 * t40 * t10
      t61 = t7 * (t14 * (-t60 * t19 * t42 * t55 + t59 * t13 * t41 * (-t1
     #7 * t13 * t24 + 1) * t19 * t10 * t15) + t57 * (t59 * t48 * t4 ** 2
     # * t40 * t54 * t42 * t15 * t55 + t48 * (t11 * (t58 + t38) + t16 * 
     #(t2 * t38 ** 2 * t8 ** 2 * t40 * t42 + t24 * t39)) * t54 * t15)) +
     # t47 * t43 * t1 ** 2 * t2 * t34 * t36
      t62 = t2 * t22
      t63 = t4 * t13
      t64 = t63 - t62
      t23 = t23 * t39
      t24 = t38 * t24
      t39 = t63 * t38
      t32 = t48 * t1 * ((t15 * (t23 * t64 * t41 + t24 * (-t63 + t62)) + 
     #t38 + t58 + t39 * t64 * t41 * t55) * t19 * t14 + t43 * t2 * (-t2 *
     # (t27 * t32 + t52) * t34 * t26 - t44))
      t56 = t6 * t56
      t4 = t48 * t20 * t40 * t15 * t43 * (t56 * t31 * t35 + t16 * t38 * 
     #t41 * (-t2 * t7 + t25 * t4))
      t3 = t3 * t2 * ((t55 * (-t63 * t53 * t16 * t54 * t41 + t60 * t8 * 
     #t33 * t54 * t42) + t10 * t33 * (-t59 * t8 * t41 + t56) * t54 * t13
     # * t15) * t57 * t7 + t43 * t1 * t25 * t20 * t36 * (-t16 * t34 * t3
     #3 + t10))
      t6 = t51 * t19 * t14 * t7 * (t39 * t55 * t10 * t41 - t50 * (t13 * 
     #t6 * t11 + t2) + (t23 * t41 - t24) * t10 * t15)
      t8 = 4
      ret = 16 * t2 * (-t12 * t15 * t16 * t38 * t40 * t7 * t10 * t42 + t
     #48 * (t44 * t31 * t15 + t45 * t25 * (-t16 * t26 * t34 * t27 + 1) +
     # t17 * t26 * (-t16 * t25 * t34 + t31) * t29) * t20 * t43 - t21 * t
     #48 * t14 * t19 * t13 * t11 + t1 * (t10 * (t46 + t2) + t47 * (t30 +
     # t5)) * t36) + t8 * (t18 * t1 * t38 * t40 * t10 * t41 - t52 * t48 
     #* t38 * t31 * t26 - t53 * (t45 + t51) + t51 * (t11 * (-t1 * t35 + 
     #t2 * t33) + t37 * (t15 * t2 * t10 - t1)) * t7 + t28 * t10 * (t15 *
     # (-t49 * t22 * t41 - t37 * t47 * t31) + t50 * (t46 + t2))) - 64 * 
     #t61 - 8 * t32 - 48 * t28 * t48 * t43 * t25 * t20 * (t37 * t16 + t1
     #1) + 32 * t4 - 128 * t3 - 24 * t6 + 12 * t49 * t15 * t9 * t13 * t7
     # * t10 * t41

      hjetmass_qqbgg_triangle_pmpm_s124_mhsq_0 = ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     &  hjetmass_qqbgg_triangle_pmpm_s124_mhsq_0_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          real*16 cg

      t1 = za(i2, i4)
      t2 = zb(i3, i1)
      t3 = za(i1, i2)
      t4 = za(i2, i3)
      t5 = zb(i2, i1)
      t6 = za(i1, i4)
      t7 = zb(i4, i1)
      t8 = zb(i4, i2)
      t9 = t3 * t5
      t10 = t6 * t7
      t11 = t1 * t8 + t10 + t9
      t12 = za(i1, i3)
      t13 = zb(i3, i2)
      t14 = za(i3, i4)
      t15 = zb(i4, i3)
      t16 = t12 * t2
      t17 = t4 * t13
      t18 = t14 * t15
      t19 = t17 + t18 + t16
      if ( real(t19) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t17 = cg * sqrt(t19 ** 2) + t16 + t17 + t18
      t17 = 0.1q1 / t17
      t18 = -0.2q1 * t4 * t2 * t11 * t17 + t1 * t7
      t9 = -0.2q1 * t16 * t11 * t17 + t10 + t9
      t10 = t12 * (t13 * t18 + t2 * t9)
      t16 = t1 * t13 + t2 * t6
      t19 = 0.2q1 * t12
      t20 = t2 * (t12 * t9 + t4 * (-t19 * t13 * t11 * t17 + t6 * t8))
      t21 = 0.2q1 * t14 * t2 * t11 * t17 - t1 * t5
      t22 = t1 * t15 + t2 * t3
      t8 = t19 * t11 * t15 * t17 - t3 * t8
      t5 = 0.1q1 / t5
      t19 = 0.1q1 / t4
      t10 = 0.1q1 / t10
      t21 = 0.1q1 / t21
      t23 = 0.1q1 / t3
      t7 = 0.1q1 / t7
      t24 = 0.1q1 / t15
      t8 = 0.1q1 / t8
      t20 = 0.1q1 / t20
      t25 = 0.1q1 / t14
      t26 = t2 * t20
      t10 = t12 * t10
      t12 = t9 * (t25 * t8 + t26)
      t27 = t7 * t23 * t5
      t28 = t27 * t2
      t11 = 0.1q1 / t11
      ret = 0.8q1 * t28 * t1 * (t12 + t2 * (t21 * t24 + t10) * t19 * t18
     #) - 0.48q2 * t28 * t9 * t22 * t17 * (t26 * t14 + t8) + 0.16q2 * t2
     #8 * t17 * (t2 * (t22 * t24 * (-t14 * t18 * t19 * t21 + 0.1q1) + t1
     #0 * t18 * (-t14 * t22 * t19 + t16)) + t12 * t16 * t4) + 0.128q3 * 
     #t2 * t1 * t22 * t17 * t23 * t11 * (t14 * t19 * t5 - t7) - 0.32q2 *
     # t27 * t17 * t9 ** 2 * t4 * (-t16 * t15 * t25 * t8 ** 2 + t14 * t2
     # ** 2 * t20 ** 2 * (-t13 * t22 + t2 * (t13 * t3 - t15 * t6))) - 0.
     #64q2 * t1 ** 2 * t2 * t23 * t19 * t5 * t11

      hjetmass_qqbgg_triangle_pmpm_s124_mhsq_0_rat = ret/32q0*(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_triangle_pmpm_s12_s123_0
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt
          real*16 cg

      t1 = za(i2, i4)
      t2 = zb(i3, i2)
      t3 = za(i3, i4)
      t4 = za(i1, i2)
      t5 = zb(i4, i2)
      t6 = za(i1, i3)
      t7 = zb(i2, i1)
      t8 = zb(i4, i3)
      t9 = zb(i3, i1)
      t10 = za(i2, i3)
      t11 = t6 * t9
      t12 = t10 * t2
      t13 = t12 + t11
      if ( real(t13) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t13 = cg * sqrt(t13 ** 2) + t11 + t12
      t14 = 0.1q1 / t13
      t15 = 2 * t6
      t16 = t15 * t7
      t17 = t4 * (t16 * t8 * t14 - t5)
      t18 = t4 * t7
      t15 = t18 * (-t15 * t9 * t14 + 1)
      t16 = t9 * (t15 * t6 - t17 * t3) - t12 * t16 * t4 * t9 * t14
      t19 = za(i1, i4)
      t20 = t1 * t2
      t21 = t19 * t9
      t22 = zb(i4, i1)
      t17 = 0.1q1 / t17
      t16 = 0.1q1 / t16
      t23 = 0.1q1 / t8
      t13 = 0.1q1 / t13
      t24 = t12 + t11
      t25 = t16 ** 2
      t26 = t16 * t25
      t27 = t17 ** 2
      t28 = t3 ** 2
      t29 = t2 ** 2
      t30 = t4 ** 2
      t31 = t4 * t30
      t32 = t9 ** 2
      t33 = t9 * t32
      t34 = t19 * t22
      t35 = t1 * t5
      t36 = t18 * t13
      t37 = t23 * t15
      t38 = t37 * t33
      t39 = t36 * t33
      t40 = t33 * t15
      t41 = 0.1q1 / t3
      t42 = 0.1q1 / t7
      t43 = t34 + t35
      t44 = mt ** 2
      t45 = t3 * t7
      t46 = t3 * t15 * t16
      t47 = t32 * t16
      t48 = t41 * t17
      t49 = t12 * t23 + t3
      t50 = t1 ** 2
      t51 = t10 * t19
      t52 = t1 * t6
      t53 = t15 * t42
      t54 = t11 * t23
      t55 = t9 * t3
      t56 = t9 * t16
      t57 = t1 * t9
      t5 = t2 * (t45 * t38 * t12 * t31 * t26 + t47 * t13 * (t1 * t29 * t
     #7 * t10 ** 2 * t16 * t23 + t12 * t16 * (-t55 * t37 + (t54 + t3) * 
     #t7 * t1) + t55 * ((-t54 - t3) * t16 * t15 + t23)) * t30 + t57 * t4
     #2 * t23 * (t56 + t48) * t44 + t10 * (t32 * (t20 * t15 * t13 * t49 
     #* t25 - t20 * t23 * t13 * t16) + t33 * (-t19 * t23 * t13 * t16 + t
     #13 * t15 * (t19 * t3 + (t52 + t51) * t23 * t2) * t25 + t37 * t2 * 
     #t3 * t42 * (t19 ** 2 * t22 ** 2 + t5 ** 2 * t50) * t26) + t37 * t6
     # * t32 ** 2 * t19 * t25 * t13 + t53 * t2 * t8 * t17 * t27) * t4)
      t7 = t53 + t4
      t19 = t43 * t42
      t1 = t47 * t23 * t2 * (t16 * (t21 * t15 * (t19 + t4) + t20 * (t4 *
     # (t18 + t15) + t34 * t7 + t35 * t7)) * t10 + t46 * t9 * t4 * (-t19
     # - t4) + t11 * t4 * t1 * t14)
      t7 = t2 * (t17 * (t4 * t32 * t42 * t23 + t41 * t50) + t56 * (t3 * 
     #t4 * t32 * t42 * t23 + t50 + t57 * (t41 * (t52 - t51) + t53) * t23
     #))
      ret = -64 * t13 * t10 * t30 * t29 * (t27 * (-t36 * t9 + (-t36 - 1)
     # * t17 * t15 * t8) + t28 * (t39 * t25 + t40 * (t18 + t35 + t34) * 
     #t26) + t3 * (t39 * t23 * t24 * t25 + t38 * (t34 * t24 + t35 * t24 
     #+ t18 * (-t6 ** 2 * t32 * t13 + t11 + t12 * (-t12 * t13 + 1))) * t
     #26) - t36 * t3 * t28 * t33 * t15 * t8 * t26) - 32 * t12 * (t2 * (-
     #t35 * t34 * t38 * t3 * t26 * t42 * t4 + t9 * (-t13 * t27 + (-t13 *
     # t43 * t23 * t25 - t37 * t43 * t26) * t32 * t3) * t30 - t45 * t33 
     #* t25 * t23 * t13 * t31) + t42 * t14 * (t20 + t21) * t44 * (t15 * 
     #t41 * t27 + t47 * (-t46 + t23) + t48 * t9 * t23)) + 128 * t45 * t4
     #0 * t13 ** 2 * t26 * t10 * t31 * t29 * (t11 * t49 + t12 * t3) + 16
     # * t5 - 8 * t1 - 4 * t7

      hjetmass_qqbgg_triangle_pmpm_s12_s123_0 = ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_triangle_pmpm_s12_s124_0
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt
          real*16 cg

      t1 = za(i1, i2)
      t2 = zb(i2, i1)
      t3 = za(i1, i4)
      t4 = zb(i4, i1)
      t5 = za(i2, i4)
      t6 = zb(i4, i2)
      t7 = t3 * t4
      t8 = t5 * t6
      t9 = t7 + t8
      if ( real(t9) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t9 = cg * sqrt(t9 ** 2) + t7 + t8
      t10 = zb(i4, i3)
      t11 = za(i3, i4)
      t12 = 0.1q1 / t9
      t13 = -2 * t7 * t12
      t14 = t1 * t2
      t15 = t14 * (1 + t13)
      t16 = za(i2, i3)
      t17 = t2 * (-2 * t1 * t11 * t4 * t12 - t16)
      t13 = t8 * t14 * t13 + t3 * (t10 * t17 + t15 * t4)
      t18 = zb(i3, i1)
      t19 = zb(i3, i2)
      t20 = za(i1, i3)
      t21 = t18 * t3 + t19 * t5
      t22 = 0.1q1 / t11
      t17 = 0.1q1 / t17
      t9 = 0.1q1 / t9
      t13 = 0.1q1 / t13
      t23 = t16 * t19
      t24 = t18 * t20
      t25 = t14 + t23 + t24
      t26 = t13 ** 2
      t27 = t13 * t26
      t28 = t4 ** 2
      t29 = t17 ** 2
      t30 = t10 ** 2
      t31 = t5 ** 2
      t32 = t5 * t31
      t33 = t2 ** 2
      t34 = t3 ** 2
      t35 = t34 ** 2
      t36 = t3 * t34
      t37 = t8 * t22
      t38 = t14 * t9
      t39 = t25 * t22
      t40 = t9 * t12
      t41 = 0.1q1 / t10
      t42 = t24 + t23
      t43 = t9 ** 2
      t44 = mt ** 2
      t45 = t1 * t33
      t46 = t45 * t5 * t43
      t47 = t46 * t26
      t48 = t12 * t22
      t49 = t48 * t5 * t4
      t50 = t44 * t4 * t21 * t12 ** 2
      t51 = t31 * t4
      t52 = t51 * t2
      t53 = 0.1q1 / t1
      t54 = t2 * t10
      t55 = t4 * t19
      t56 = t55 - t54
      t57 = t18 ** 2
      t58 = t1 ** 2
      t59 = t6 * t18
      t60 = t10 * t18
      t61 = t45 * t40
      t21 = t48 * t44 * t21 * t53
      t48 = t9 * t5
      t62 = t3 * t13
      t16 = t51 * (t17 * (-t21 * t41 + (t11 * t4 * t12 * t29 + t17 * t9)
     # * t5 * t33) + t36 * (t49 * t10 * t33 * (-t16 ** 2 * t19 ** 2 - t5
     #7 * t20 ** 2 - t33 * t58) * t27 + t61 * t4 * (t60 + (t59 - t54 + t
     #55) * t22 * t5) * t26) + t48 * t33 * (t10 * (t1 * t56 * t12 + t39)
     # + t37 * t1 * t12 * t56) * t26 * t34 + t62 * t21 + t61 * t35 * t18
     # * t28 * t22 * t26)
      t19 = -t55 + t54
      t21 = t4 * t2
      t51 = t62 * t22
      t54 = t17 * t41
      t61 = t22 * t18
      t19 = t31 * (t34 * (t21 * t18 * t22 * t9 * t13 + t21 * (t60 * t38 
     #+ (t18 * (t12 * t19 * t20 + t38 * t6) + t14 * t12 * t19 + t23 * t1
     #2 * t19) * t22 * t5) * t26) + t51 * (-t44 * t18 * t53 + t48 * t2 *
     # t56) + t54 * t44 * t18 * t53 * t22 + t61 * t28 * t2 * (-t12 * t42
     # + t14 * (-t12 + t9)) * t26 * t36)
      ret = 128 * t40 * t32 * t28 * t2 * t33 * t1 * (t36 * (-t14 * t11 *
     # t9 * t27 * t10 * t30 + t37 * (t14 * (-t8 * t9 + 1) + t23 + t24) *
     # t27 * t10 + t25 * t27 * t30) + t39 * t4 * t27 * t10 * t35 + t11 *
     # t17 * t29 * (t38 + 1) - t38 * t3 * t35 * t28 * t10 * t22 * t27) -
     # 64 * t52 * (t29 * (t50 * t41 - t46) + t34 * (t47 * t30 + (t45 * t
     #31 * t6 * t22 * t43 - t50) * t26 * t10) + (t47 * t4 * t22 + t49 * 
     #t2 * (t14 * t42 + t24 * t23) * t27) * t10 * t36) + 32 * t16 - 256 
     #* t27 * t43 * t12 * t10 * t32 * t28 * t36 * t33 ** 2 * t58 * (t8 *
     # t10 + t7 * (t37 + t10)) + 16 * t19 + 8 * t61 * t52 * t13 * t34 * 
     #(-t13 * t25 + t12) - 4 * t5 * (t57 * (t54 - t62) + t51 * t18 * (-t
     #15 * t53 + t41 * (-t59 + t55)) * t5 + t2 * t53 * (-t62 * t10 + t17
     #) * t22 * t31)

      hjetmass_qqbgg_triangle_pmpm_s12_s124_0 = -ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_triangle_pmpm_s34_0_0
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt

      t1 = za(i1, i3)
      t2 = za(i2, i4)
      t3 = zb(i4, i1)
      t4 = za(i1, i4)
      t5 = zb(i4, i2)
      t6 = t2 * t5
      t7 = t4 * t3
      t8 = t7 + t6
      t9 = za(i2, i3)
      t10 = zb(i3, i1)
      t11 = zb(i3, i2)
      t12 = t10 * t5
      t13 = t11 * t3
      t14 = t13 - t12
      t15 = t3 ** 2
      t16 = t9 ** 2
      t17 = za(i1, i2)
      t18 = zb(i2, i1)
      t19 = t1 * t3
      t5 = t5 * t9 + t19
      t5 = 0.1q1 / t5
      t20 = 0.1q1 / t18
      t21 = 0.1q1 / t17
      t22 = t5 ** 2
      ret = 16 * za(i3, i4) * (t1 * (t16 * t10 * t14 + t2 * t15 * t8) - 
     #t9 * (-t16 * t11 * t14 + t4 * t15 * t8) - t17 * t18 * (t9 * (t3 * 
     #(2 * t7 + t6) - t9 * (2 * t13 - t12)) - t19 * (t10 * t9 + t2 * t3)
     #)) * zb(i4, i3) * t21 * t20 * t5 * t22

      hjetmass_qqbgg_triangle_pmpm_s34_0_0 = ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_triangle_pmpm_s34_mhsq_s12
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt

      t1 = za(i1, i2)
      t2 = zb(i2, i1)
      t3 = za(i1, i3)
      t4 = zb(i3, i1)
      t5 = za(i2, i3)
      t6 = zb(i3, i2)
      t7 = za(i1, i4)
      t8 = zb(i4, i1)
      t9 = za(i2, i4)
      t10 = zb(i4, i2)
      t11 = t3 * t4
      t12 = t5 * t6
      t13 = t7 * t8
      t14 = t9 * t10
      t15 = t11 + t12 + t13 + t14
      t16 = za(i3, i4)
      t17 = zb(i4, i3)
      t15 = -4 * t1 * t16 * t2 * t17 + t15 ** 2
      t15 = t11 - sqrt(t15) + t12 + t13 + t14
      t18 = 2 * t1 * t2
      t19 = t15 + t18
      t20 = (0.1q1 / 0.2q1)
      t18 = t20 * t15 ** 2 - t18 * t16 * t17
      t21 = 2 * t16 * t17 + t15
      t22 = (0.1q1 / 0.4q1)
      t23 = t15 ** 2
      t24 = t1 * t16
      t25 = t24 * t2 * t17
      t26 = t23 * t22 - t25
      t26 = 0.1q1 / t26
      t27 = t15 * t22
      t24 = t24 * t20
      t28 = t24 * t4 - t27 * t9
      t29 = t15 * t2 * t26
      t30 = t29 * t28
      t31 = -t20 * t9 * t2 * t17 + t27 * t4
      t32 = t15 * t1 * t26
      t33 = t32 * t31
      t34 = t11 + t13
      t23 = t23 * t26
      t25 = t25 * t20 * t15 * t26
      t35 = t22 * t23 * t34 - t25
      t36 = t4 * t7 + t6 * t9
      t37 = t15 * t16
      t38 = t37 * t17 * t26
      t39 = t38 * t36
      t40 = t4 * t5 + t8 * t9
      t41 = t29 * t1
      t42 = t41 * t40
      t43 = t10 * t7 + t3 * t6
      t44 = t23 * t22
      t45 = t44 * t1 * t2
      t34 = -t20 * t41 * t34 + t45
      t46 = t15 * t17 * t26
      t47 = t34 * t35
      t48 = -t42 * t23 * t43 / 8
      t49 = t48 + t30 * t46 * (t24 * t10 - t27 * t3) + t47
      t50 = t24 * t6 + t27 * t7
      t51 = t46 * t50
      t36 = t23 * t36
      t26 = t37 * t26
      t31 = -t26 * t31
      t28 = -t46 * t28
      t37 = t23 * t40
      t6 = t20 * t7 * t2 * t17 + t27 * t6
      t7 = -t32 * t6
      t40 = -t37 * t41 * t43 / 8
      t43 = t7 * t26 * (t20 * t5 * t2 * t17 + t27 * t8) + t47 + t40
      t3 = t40 + t31 * t32 * (t20 * t3 * t2 * t17 - t27 * t10) + t47
      t10 = t12 + t14
      t32 = -t20 * t41 * t10 + t45
      t13 = t13 + t14
      t14 = t22 * t23 * t13 - t25
      t5 = t48 - t51 * t29 * (t24 * t8 + t27 * t5) + t47
      t8 = t11 + t12
      t11 = t22 * t23 * t8 - t25
      t10 = t22 * t23 * t10 - t25
      t12 = t44 * t16 * t17
      t13 = -t20 * t38 * t13 + t12
      t8 = -t20 * t38 * t8 + t12
      t12 = t29 * t50
      t3 = 0.1q1 / t3
      t18 = 0.1q1 / t18
      t22 = 0.1q1 / t2
      t15 = 0.1q1 / t15
      t5 = 0.1q1 / t5
      t23 = 0.1q1 / t17
      t24 = 0.1q1 / t43
      t25 = 0.1q1 / t16
      t27 = 0.1q1 / t49
      t29 = 0.1q1 / t1
      t38 = -t34 - t32 - t11
      t40 = t34 + t32 + t14
      t41 = -t35 - t10 - t13
      t43 = t35 + t10 + t8
      t44 = t32 + t14
      t45 = t30 ** 2
      t46 = t27 ** 2
      t48 = t27 * t46
      t49 = t24 ** 2
      t50 = t24 * t49
      t52 = t18 ** 2
      t53 = t51 ** 2
      t54 = t21 ** 2
      t55 = mt ** 2
      t56 = t5 ** 2
      t57 = t5 * t56
      t58 = t3 ** 2
      t59 = t3 * t58
      t60 = t34 * t48
      t61 = t40 * t48
      t62 = t10 * t21
      t63 = t55 * t22
      t64 = t63 * t15
      t65 = t64 * t23 * t25
      t66 = t65 * t46 * t29
      t16 = t16 * t17
      t17 = t16 * t19 ** 2
      t67 = t8 * t56
      t68 = t17 * t33 * t22
      t69 = t42 * t35
      t70 = t41 * t58
      t71 = t43 * t49
      t72 = t7 * t34
      t73 = t72 * t31
      t74 = t4 * t9
      t75 = t39 * t34
      t76 = t75 * t33
      t77 = t76 * t3
      t63 = t63 * t29
      t78 = t19 * t33
      t79 = t38 * t56
      t80 = -t40 * t35 * t46 - t79 * t35 - t5
      t81 = -t10 - t8
      t82 = t32 + t11
      t83 = t35 ** 2
      t84 = t35 * t83
      t85 = t31 ** 2
      t86 = t34 ** 2
      t87 = t81 * t57
      t88 = t28 * t56
      t89 = t42 * t53
      t90 = t72 * t85
      t91 = t83 * t28
      t92 = t13 * t48
      t93 = t51 * t30
      t94 = t93 * t31
      t95 = t7 * t31
      t96 = t86 * t31
      t97 = t96 * t58
      t98 = t30 * t24
      t99 = t82 * t49
      t100 = t31 * t21
      t101 = t78 * t16
      t102 = t19 * t52
      t76 = t102 * (t101 * (t80 * t51 * t45 + t99 * t90) * t29 * t22 + t
     #100 * (-t76 * t24 + (t98 + t97) * t7 * t28))
      t103 = -t35 - t10
      t104 = t103 * t49 - t70
      t105 = t69 * t51
      t106 = t1 * t2
      t107 = t106 * t23 * t25
      t108 = t44 * t58
      t109 = t93 * t35
      t110 = t109 * (t107 * (-t28 * t13 * t46 + (-t35 * t48 + t57 * t8) 
     #* t10 * t51 * t42) * t52 * t30 * t54 - t66 * t12 * t42 * t28 + (t1
     #07 * t46 * (t105 * t27 + t28) * t18 + t102 * (t46 * (t13 * t33 + t
     #14 * t28) - t67 * t33 + t61 * t105)) * t30 * t21)
      t101 = t102 * t90 * (t21 * (t104 * t33 + t28 * (-t99 + t108)) + t1
     #01 * (t34 * t49 - t108) * t29 * t22)
      t108 = -t56 + t46
      t111 = t34 + t32
      t112 = t3 - t24
      t113 = t5 - t27
      t114 = t108 * t83
      t115 = t30 * t3
      t116 = t96 * t49
      t117 = t35 * t10
      t118 = t69 * t19 * t45 * t53 * t57 * t111 + t19 * (t45 * (t114 * t
     #33 * t51 + t69 * (t56 * (t5 * (-t35 * t82 + t38 * t8) + 1) - t46) 
     #* t53) + t77 * t31 - t95 * (t115 + t116) * t28 + t94 * t33 * t113)
     # * t18 * t21 + t117 * t107 * t45 * t51 * (-t92 * t42 * t51 + t88) 
     #* t18 * t54
      t119 = t17 * t32
      t55 = (t33 * (t31 * (-t75 * t55 * t15 * t23 * t24 * t25 + t17 * t1
     #12 * t52 * t7 * t30) - t17 * t86 * t7 * t58 * t52 * t85) + t69 * t
     #45 * t53 * t56 * (t119 * t11 * t52 * t5 + t55 * t15 * t23 * t25)) 
     #* t29 * t22
      t1 = t18 * t118 + t31 * (-t78 * t73 * t52 * t21 * t8 * t49 + (t63 
     #* (-t74 * t7 * t24 + t77 * t15) + (-t73 * t28 * t18 * t49 * t21 + 
     #t72 * (t71 + t70) * t52 * t28 * t31 * t54) * t2 * t1) * t25 * t23)
     # + t45 * (t69 * (t19 * (t18 * (t11 * t57 - t61) + t62 * (t38 * t57
     # + t60) * t52) - t66 - t60 * t17 * t22 * t29 * t44 * t52) * t53 + 
     #t52 * ((-t10 * t46 + t67) * t25 * t23 * t28 * t35 * t54 * t2 * t1 
     #+ t68 * t27 * t29) * t51) + (t63 * t74 * (t95 * t3 - t93 * t5) + (
     #(t45 * (-t91 * t46 * t51 - t92 * t89 * t83) + t94 * t28 * t27 + t8
     #5 * t28 * t7 * t3) * t52 * t54 + (t45 * (t35 * (-t88 * t51 + t42 *
     # (t10 * t48 + t87) * t53) - t89 * t57 * t83) + t90 * t28 * t58) * 
     #t18 * t21) * t2 * t1) * t25 * t23 + t76 + t110 + t101 + t55
      t2 = -t58 + t49
      t38 = t107 * t21
      t55 = t18 * t19
      t60 = t31 * t5
      t61 = t107 * t54
      t76 = t61 * t18
      t77 = t34 + t32
      t89 = t32 ** 2
      t90 = t10 ** 2
      t94 = t8 ** 2
      t101 = t11 ** 2
      t110 = t7 ** 2
      t118 = t13 ** 2
      t120 = t30 * t29
      t6 = -t33 * t26 * t6
      t26 = t6 * t58
      t121 = t59 - t50
      t9 = t9 ** 2
      t122 = t93 * t69
      t123 = t33 * t34 * t36
      t124 = t17 * t29 * t22
      t125 = t124 * t31
      t126 = t102 * t21
      t127 = t61 * t52
      t128 = t127 * (t42 * t39 * t51 * t27 * (-t30 * t83 * t27 + t31) + 
     #t85 * ((-t90 - t83) * t50 * t34 + t70) * t110 * t37)
      t4 = t4 ** 2
      t129 = t14 ** 2
      t130 = t36 * t34 * t31
      t131 = t34 * t31
      t132 = t9 * t29 * t25
      t133 = t95 * (-t112 * t23 * t22 * t4 - t132 * t3)
      t4 = t93 * (t113 * t23 * t22 * t4 + t132 * t113)
      t132 = t107 * t34 * t37 * t121 * t110 * t85
      t134 = t126 * (-t130 * t28 * t3 + (-t30 * t8 - t31 * t82) * t31 * 
     #t49 * t110 * t37 + t122 * (t36 * t8 - t75) * t56)
      t111 = t125 * t52 * t3 * (t123 + (-t115 * t111 + t131 * (t129 + t8
     #9) * t58) * t110 * t37)
      t4 = t18 * (t37 * (t55 * t73 * t21 * (t39 * t11 * t49 + t26) + t31
     # * (t107 * t100 * ((t35 * t49 + (t118 + t83 + t90) * t34 * t59 - t
     #34 * t94 * t50) * t18 * t21 - t49 + t58) + t55 * (t100 * t58 * (t3
     #2 + t14) + t16 * (t120 * t49 - t31 * (t101 + t89) * t29 * t50) * t
     #22 * t34 * t19)) * t110) + (t68 * t36 * t18 * t29 * (-t5 + t27) + 
     #t76 * t56 * (t10 + t8) * t51 * t42 * t39 + (t55 * t46 * t77 - t107
     # * t56) * t51 * t42 * t39 * t21) * t35 * t30 + t55 * t100 * t36 * 
     #t34 * t28 * t24) + t126 * (-t122 * t39 * t56 * t82 + (t30 * t104 *
     # t110 - t6 * t72 * t49) * t37 * t31) + t128 + t95 * t9 * t29 * t25
     # * t24 + t125 * (-t123 * t24 + (t30 * (-t14 * t58 + t99) + t121 * 
     #t31 * t34 * t86) * t110 * t37) * t52 + t111 + t132 + t134 - t127 *
     # t85 * t37 * t110 * t49 * t81 + t133 + t4
      t9 = t10 + t13
      t68 = t46 * t9
      t72 = t30 * t35
      t99 = t72 * t39
      t4 = t4 + (t95 * t64 * t75 * t37 * t11 * t49 * t29 + (t120 * t67 *
     # t64 * t35 * t36 + t106 * ((t30 * (-t68 * t35 + t56 * t83) - t60) 
     #* t52 * t39 * t54 + t99 * t18 * t46 * t21)) * t51 * t42) * t25 * t
     #23 + t55 * t30 * t37 * t110 * t31 * t2
      t67 = t127 * t31 * t39
      t83 = t21 * t13
      t104 = t65 * t35
      t111 = t36 * t42
      t122 = t120 * t55 * t16
      t123 = t7 * t37
      t67 = (t21 * (-t107 * t75 * t31 * t2 * t18 + t102 * (-t116 * t39 +
     # t98 * t39 + t131 * t58 * (-t13 * t36 + t39 * t44))) + t65 * t36 *
     # t31 * t29 * t112 + t67 * (t70 * t34 + t71 * t34 + t3)) * t7 * t37
     # + t111 * (-t104 * t39 * t29 * t113 + t124 * t93 * t52 * t27 + t93
     # * (t35 * t83 * t46 * t52 + t35 * t56 * t18) * t19) + t123 * (-t13
     #0 * t65 * t13 * t29 * t58 - t67 * t24 + t126 * (t39 * (-t131 * t32
     # * t49 - t115 + t97) + t31 * t112 * t36)) + t111 * t55 * ((t113 * 
     #t31 + t30 * (t117 * t108 + t114)) * t18 * t21 + t122 * t80 * t22 -
     # t72 * t46) * t51
      t17 = t120 * t17 * t22
      t65 = t65 * t29
      t70 = t126 * t27
      t43 = t43 * t50
      t71 = t64 * t49 * t29
      t80 = t16 * t19 * t22
      t97 = t102 * t85
      t98 = t18 * t48
      t111 = t52 * t21
      t26 = -t65 * t26 * t95 * t37 + t111 * (-t19 * t58 + t19 * (-t10 * 
     #t11 - t8 * t82) * t50 - t107 * t62 * t13 * t59) * t110 * t37 * t85
      t83 = t97 * t59 * (-t80 * t32 * t29 + t83) * t110 * t37 * t86
      t47 = t30 * (t107 * (t41 * t46 * t52 * t31 * t54 + t72 * (t48 - t5
     #7)) + t55 * (t30 * t56 + (t31 * (t34 * t46 + t79) + t68 * t30) * t
     #18 * t21 + t122 * t34 * t56 * (-t47 * t5 + 1) * t22)) * t53 + t126
     # * t109 * t12 * t28 * t46
      t68 = t36 * t28
      t79 = t110 * t37
      t62 = t34 * (-t119 * t37 * t110 * t14 * t59 * t22 * t29 * t85 - t6
     #1 * t39 * t28 * t24 * t31) + t79 * t19 * (t59 * (-t80 * t14 * t29 
     #+ t62) - t43 * t21) * t85 * t86 + t72 * t19 * t21 * (t33 * t39 * t
     #5 + t68 * t27)
      t17 = t42 * t47 + t52 * t62 + t34 * t26 + t42 * (t93 * t39 * (-t65
     # * t5 + t70) + t18 * t30 * (t38 * (t108 * t31 + (t98 * t118 * t35 
     #+ t98 * t84) * t30 * t21) + t17 * t35 * t18 * (t48 * (t129 + t89 +
     # t86) + t57 * (-t101 - t89)) + t46 * t19 * (t100 * t18 * t44 - t30
     #)) * t53) + t30 * (t42 * (t39 * t27 * (-t70 * t35 * t14 + t65) * t
     #51 + t52 * (t17 * t56 * t82 + t61 * (-t103 * t56 * t31 + t72 * t90
     # * t48)) * t53) - t104 * t36 * t28 * t29 * t113) + t85 * ((t106 * 
     #(-t41 * t59 - t43) * t18 * t21 + t71) * t25 * t23 + t55 * (t11 * t
     #50 + (t50 * (t103 * t32 - t11 * t35) + t35 * t32 * t59) * t18 * t2
     #1)) * t110 * t37 * t34 + t97 * (t80 * t50 * t29 * t82 + t21 * t35 
     #* t59) * t110 * t37 * t86 + t99 * t52 * t21 * t27 * (t38 * t28 - t
     #78) + t83
      t26 = t79 * t85
      t6 = t42 * (-t109 * t66 * t14 * t39 + t52 * t45 * (-t124 * t46 * t
     #40 - t61 * t84 * t57) * t53) + t17 + t26 * ((-t64 * t58 * t29 + t1
     #06 * (t50 * (-t103 * t8 + t117) - t9 * t35 * t59) * t52 * t54) * t
     #25 * t23 + t55 * (t59 * ((-t14 * t41 + t32 * t9) * t18 * t21 - t32
     # - t14) + t32 * (t55 * t16 * t11 * t22 * t29 + 1) * t50 + t21 * t1
     #8 * t49)) * t34 + t111 * t5 * t30 * (t42 * (t5 * (t19 * t30 * t81 
     #+ t38 * (t31 * t8 + t72 * (-t90 - t94) * t5)) * t53 - t19 * (t28 *
     # t35 * t12 * t5 + t39) * t51) - t68 * t19 * t35) - t26 * t55 * t12
     #1 * t86 - t99 * t127 * t28 * t5 + t131 * t23 * t25 * (t106 * t54 *
     # t39 * t28 * t3 * t52 + t71 * t123 * t6)
      t9 = -t58 + t49
      t3 = t3 - t24
      t17 = t73 * t37
      t3 = t36 * (t19 * (t17 * t9 * t18 + t111 * (t17 * (t49 * (-t35 - t
     #10 - t8) + t58 * (t35 + t10)) + t69 * t39 * (t5 - t27))) + t65 * t
     #75 * t37 * t3 + t123 * t124 * (t30 * t3 + t131 * (t49 * (t32 + t11
     #) + t58 * (-t32 - t14)) + t96 * t9) * t52)
      t5 = 16
      ret = t20 * t126 * t75 * t36 * t37 * t112 + t5 * (t102 * t51 * t35
     # * t45 * t21 * (-t105 * t34 * t57 + (-t56 + t46) * t10 * t33 - t88
     # * (t34 + t32 + t11) + t28 * t77 * t46) + t45 * (t69 * t18 * (t38 
     #* (-t87 * t21 * t35 * t18 + t92) + t55 * (t48 * t21 * (t10 * t44 +
     # t13 * t40) + t16 * (-t32 * t14 * t48 + t82 * t57 * t34) * t29 * t
     #22 * t19)) * t53 + t18 * (t107 * t91 * t54 * t18 * t56 + t55 * t11
     #3 * t28 * t21 - t78 * t35 * t108) * t51) + t1 + t18 * t7 * t85 * (
     #t78 * (-t112 * t18 * t21 + t2 * t34) - t76 * t28 * t24) + t93 * t2
     #3 * t25 * (-t60 * t106 * t54 * t28 * t52 + t63 * (t88 * t69 * t12 
     #* t15 + t74 * t27))) - 8 * t6 + 4 * t4 - 2 * t67 + t3

      hjetmass_qqbgg_triangle_pmpm_s34_mhsq_s12 = ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     &  hjetmass_qqbgg_triangle_pmpm_s34_mhsq_s12_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          real*16 cg

          cg = 1q0

      t1 = za(i2, i4)
      t2 = zb(i3, i1)
      t3 = za(i1, i3)
      t4 = za(i2, i3)
      t5 = zb(i3, i2)
      t6 = za(i1, i4)
      t7 = zb(i4, i1)
      t8 = zb(i4, i2)
      t9 = t3 * t2
      t10 = t4 * t5
      t11 = t6 * t7
      t12 = t1 * t8
      t13 = t9 + t10 + t11 + t12
      t14 = za(i1, i2)
      t15 = za(i3, i4)
      t16 = zb(i2, i1)
      t17 = zb(i4, i3)
      t18 = t14 * t16
      t13 = -0.4q1 * t18 * t15 * t17 + t13 ** 2
      t13 = cg * sqrt(t13) + t10 + t11 + t12 + t9
      t19 = 0.1q1 / 0.2q1
      t20 = t19 * t13
      t21 = t6 * t16 * t17 + t20 * t5
      t22 = -t14 * t21
      t23 = t1 * t16 * t17 - t20 * t2
      t24 = t15 * t23
      t25 = t11 + t9
      t26 = t13 ** 2
      t27 = t25 * t16
      t28 = t15 * t5
      t29 = t26 * t14
      t30 = t4 * t2
      t31 = t1 * t7
      t32 = t31 + t30
      t33 = t3 ** 2 * t2 ** 2
      t34 = t6 ** 2 * t7 ** 2
      t35 = t5 * t32 * t3
      t32 = t8 * t32 * t6
      t36 = t18 * t13
      t37 = t17 ** 2
      t25 = t25 * t14
      t38 = t25 * t17
      t4 = t4 * t6
      t39 = t11 * t9 * t13
      t40 = 0.1q1 / 0.4q1
      t41 = t19 * t36 * ((t11 + t18 + t10) * t17 * t15 + t33 + t34 + t35
     # + t32)
      t42 = -t40 * t29 * (-t28 * t7 + t27) - t18 * ((-t4 * t37 + t38) * 
     #t16 * t15 - t39) + t41
      t43 = t15 * t2
      t3 = t1 * t3
      t32 = t19 * t36 * ((t12 + t18 + t9) * t17 * t15 + t33 + t34 + t35 
     #+ t32)
      t29 = -t40 * t29 * (t43 * t8 + t27) - t18 * ((t3 * t37 + t38) * t1
     #6 * t15 - t39) + t32
      t33 = -t20 * t1 + t43 * t14
      t34 = t16 * t33
      t28 = t28 * t14 + t20 * t6
      t35 = t17 * t28
      t36 = t26 * t16
      t37 = t15 ** 2
      t27 = t27 * t17 * t15
      t4 = t40 * t36 * (t14 * (-t11 - t9) + t4 * t17) + t41 + t18 * (t14
     # * (t5 * t7 * t17 * t37 - t27) + t39)
      t3 = -t40 * t36 * (t3 * t17 + t25) + t32 - t18 * (t14 * (t2 * t8 *
     # t17 * t37 + t27) - t39)
      t5 = t1 * t5 + t2 * t6
      t6 = t31 + t30
      t7 = t18 * t6
      t8 = t11 + t9
      t25 = t18 * t15 * t17
      t27 = t19 * t13 * t8 - t25
      t30 = t26 * t40 - t25
      t8 = t20 * t18 - t18 * t8
      t18 = -t17 * t33
      t9 = t9 + t10
      t10 = t15 * t17
      t20 = t10 * t20
      t11 = t11 + t12
      t12 = t10 * t5
      t23 = -t14 * t23
      t31 = 0.1q1 / t42
      t3 = 0.1q1 / t3
      t29 = 0.1q1 / t29
      t17 = 0.1q1 / t17
      t14 = 0.1q1 / t14
      t32 = 0.1q1 / t16
      t4 = 0.1q1 / t4
      t33 = 0.1q1 / t15
      t36 = -t3 + t4
      t37 = t34 * t35
      t38 = t33 * t32 * t14 * t17
      t39 = -t31 + t29
      t30 = 0.1q1 / t30 ** 2
      t40 = t27 * (t3 ** 2 - t4 ** 2)
      t41 = t31 ** 2
      t42 = -t29 ** 2 + t41
      t14 = t32 * t14 * t30
      ret = 0.16q2 * t38 * t2 * t1 * (t37 * t36 + (t31 - t29) * t24 * t2
     #2) - 0.8q1 * t38 * t30 * t13 * (t40 * t7 * t35 ** 2 * t34 ** 2 + t
     #8 * t23 * t24 * t12 * t39 + t37 * (t12 * (t3 * (-t27 * (t19 * t13 
     #* t11 - t25) * t3 - 0.1q1) + t4) + t40 * t18 * t16 * t28) * t7) - 
     #0.2q1 * t14 * t6 * t5 * t13 * t26 * (t5 * t8 * t39 + (t29 * (-t8 *
     # (-t10 * t11 + t20) * t29 - 0.1q1) + t31) * t33 * t17 * t24 * t22)
     # - 0.4q1 * t26 * t14 * (t5 ** 2 * t7 * t27 * t36 + (t42 * t8 * t6 
     #* t24 ** 2 * t22 ** 2 + t27 * t5 * t34 * (-t18 * t3 + t4 * (-t35 *
     # t7 * (-t10 * t9 + t20) * t4 + t18)) + (-t42 * t15 * t21 * t23 - t
     #12 * (t19 * t13 * t9 - t25) * t41) * t8 * t6 * t24 * t22) * t33 * 
     #t17)

      hjetmass_qqbgg_triangle_pmpm_s34_mhsq_s12_rat = ret/32q0*(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_triangle_pmpp_s123_mhsq_0
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt
          real*16 cg

      t1 = za(i1, i2)
      t2 = zb(i2, i1)
      t3 = za(i1, i3)
      t4 = zb(i3, i1)
      t5 = za(i2, i3)
      t6 = zb(i3, i2)
      t7 = t1 * t2
      t8 = t3 * t4
      t9 = t5 * t6 + t7 + t8
      t10 = za(i1, i4)
      t11 = zb(i4, i1)
      t12 = za(i2, i4)
      t13 = zb(i4, i2)
      t14 = za(i3, i4)
      t15 = zb(i4, i3)
      t16 = t10 * t11
      t17 = t12 * t13
      t18 = t14 * t15
      t19 = t16 + t17 + t18
      if ( real(t19) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t18 = cg * sqrt(t19 ** 2) + t16 + t17 + t18
      t19 = t1 * t13 + t15 * t3
      t20 = 0.1q1 / t18
      t8 = -2 * t16 * t9 * t20 + t7 + t8
      t21 = -2 * t10 * t9
      t22 = t10 * t8
      t23 = t11 * (t12 * (t21 * t13 * t20 + t3 * t6) + t22)
      t6 = t21 * t15 * t20 - t1 * t6
      t21 = t5 * t4
      t24 = -2 * t12 * t9 * t11 * t20 + t21
      t25 = t1 * t11 - t15 * t5
      t26 = t10 * (t11 * t8 + t13 * t24)
      t27 = -2 * t14 * t9 * t11 * t20 - t2 * t5
      t26 = 0.1q1 / t26
      t28 = 0.1q1 / t1
      t2 = 0.1q1 / t2
      t23 = 0.1q1 / t23
      t29 = 0.1q1 / t6
      t27 = 0.1q1 / t27
      t30 = 0.1q1 / t5
      t14 = 0.1q1 / t14
      t31 = t4 * t19
      t32 = t11 * t6
      t33 = t15 * t8
      t34 = t33 - t32 - t31
      t35 = t11 * t30
      t36 = t15 * t28
      t37 = -t35 + t36
      t38 = t23 ** 2
      t39 = t23 * t38
      t40 = t11 ** 2
      t41 = t11 * t40
      t42 = t15 ** 2
      t43 = t12 * t19
      t44 = t33 * t29
      t45 = t43 * t28
      t46 = t40 * t12
      t47 = t4 * t15
      t48 = t11 * t2
      t49 = 0.1q1 / t9
      t18 = 0.1q1 / t18
      t50 = t12 ** 2
      t51 = mt ** 2
      t52 = t30 * t2
      t53 = t52 * t43 * t4
      t54 = t19 * t6
      t55 = t2 * t49
      t56 = t11 * t23
      t57 = t40 * t6
      t3 = t11 * (t9 * (t11 * (-t53 * t23 * t18 + t8 * (t4 * t13 * t2 + 
     #t15) * t30 * t18 * t38 * t19 * t50) + t40 * (t54 * t38 * t18 * t30
     # * t50 + t22 * t53 * t38 * t18)) + t55 * t47) + t52 * (t25 * (t11 
     #* (t32 * t23 - t14) + t33 * (-t14 * t29 + t56)) * t12 + t8 * (t57 
     #* t38 * (t11 * t19 - t13 * t25) + t15 * (t11 * t3 + t13 * t5) * t2
     #9 * t14 ** 2) * t50 + t11 * t25 * (t10 * t15 * t26 - t27) * t24) *
     # t20 * t28 * t51
      t5 = t32 + t31
      t5 = t43 * t5 * t30 - t5 * t6
      t52 = -t43 * t30 + t6
      t53 = t29 ** 2
      t58 = t18 * t2
      t59 = t32 * t8 * t38
      t60 = t30 * t18
      t10 = t60 * t23 * t19 * t50 * t40 * t9 * (t59 * (t17 + t16) + t58 
     #* (t10 * (t57 * t23 + t56 * t33) - t15 - t10 ** 2 * t41 * t8 * t6 
     #* t38 + t17 * t23 * (t33 + t32) - t59 * t50 * t13 ** 2) * t28 * t9
     #)
      t1 = -t60 * t44 * t48 * t9 * t12 * t14 + t1 * t4 * (t43 * t8 * t38
     # + t55) * t30 * t40 + t21 * t55 * t42 * t28 - t7 * t54 * t50 * t41
     # * t8 * t30 * t39
      t7 = -4
      ret = t7 * (-t48 * t37 * t27 * t24 + t48 * (-t46 * t6 * t30 + t32 
     #* (t36 * t12 - t4) + t47 * (t45 + t8)) * t23 + t46 * (t43 * t34 * 
     #t30 - t34 * t6) * t38 + t2 * (-t12 * t42 * t8 * t28 * t29 + t44 * 
     #(t35 * t12 - t4) + t4 * t11 * (-t43 * t29 * t30 + 1)) * t14 + t16 
     #* t15 * t24 * t2 * t37 * t26) - 32 * t3 - 8 * t58 * t12 * t9 * (t3
     #0 * (-t11 * (t31 * t29 + t11) + t31 * t33 * t53 - t42 * t8 ** 2 * 
     #t53 + t45 * t29 * t15 * (t29 * (-t33 + t31) + t11)) * t14 + t56 * 
     #t36 * (-t33 + t31 + t32) + t28 * t40 * (t16 * t5 + t17 * t5 + t33 
     #* (t16 * t52 + t17 * t52)) * t38) + 64 * t10 + 16 * t1 - 128 * t2 
     #* (t22 * t54 * t9 ** 2 * t18 ** 2 * t12 * t50 * t40 ** 2 * t13 * t
     #28 * t30 * t39 + (-t36 + t35) * t20 * t49 * t25 * t4 * t51)

      hjetmass_qqbgg_triangle_pmpp_s123_mhsq_0 = ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     &  hjetmass_qqbgg_triangle_pmpp_s123_mhsq_0_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          real*16 cg

      t1 = za(i2, i4)
      t2 = zb(i4, i1)
      t3 = za(i1, i2)
      t4 = za(i2, i3)
      t5 = zb(i4, i3)
      t6 = t2 * t3 - t4 * t5
      t7 = za(i1, i4)
      t8 = zb(i4, i2)
      t9 = za(i3, i4)
      t10 = t7 * t2
      t11 = t1 * t8
      t12 = t9 * t5
      t13 = t12 + t10 + t11
      if ( real(t13) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t11 = cg * sqrt(t13 ** 2) + t10 + t11 + t12
      t12 = zb(i2, i1)
      t13 = zb(i3, i1)
      t14 = za(i1, i3)
      t15 = zb(i3, i2)
      t16 = t3 * t12
      t17 = t14 * t13
      t18 = t15 * t4 + t16 + t17
      t11 = 0.1q1 / t11
      t19 = -0.2q1 * t1 * t18 * t2 * t11 + t13 * t4
      t20 = -0.2q1 * t9 * t18 * t2 * t11 - t12 * t4
      t10 = -0.2q1 * t10 * t18 * t11 + t16 + t17
      t16 = t7 * (t10 * t2 + t19 * t8)
      t17 = -0.2q1 * t7 * t18
      t21 = t2 * (t1 * (t17 * t8 * t11 + t14 * t15) + t10 * t7)
      t15 = t17 * t5 * t11 - t15 * t3
      t17 = 0.1q1 / t21
      t9 = 0.1q1 / t9
      t12 = 0.1q1 / t12
      t20 = 0.1q1 / t20
      t16 = 0.1q1 / t16
      t21 = 0.1q1 / t4
      t22 = 0.1q1 / t15
      t23 = 0.1q1 / t3
      t24 = t2 * t17
      t18 = 0.1q1 / t18
      ret = -0.32q2 * t23 * t11 * t21 * t12 * (t6 * (t2 * (t24 * t15 - t
     #9) + (-t22 * t9 + t24) * t10 * t5) * t1 + t2 * t6 * (t5 * t7 * t16
     # - t20) * t19 + t10 * (t5 * (t14 * t2 + t4 * t8) * t9 ** 2 * t22 +
     # t2 ** 2 * t15 * t17 ** 2 * (t2 * (t14 * t5 + t3 * t8) - t6 * t8))
     # * t1 ** 2) - 0.128q3 * t13 * t6 * t11 * t12 * t18 * (t2 * t21 - t
     #23 * t5)

      hjetmass_qqbgg_triangle_pmpp_s123_mhsq_0_rat = ret/32q0*(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_triangle_pmpp_s124_mhsq_0
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt
          real*16 cg

      t1 = zb(i4, i1)
      t2 = zb(i4, i3)
      t3 = za(i1, i2)
      t4 = zb(i3, i1)
      t5 = za(i2, i4)
      t6 = t3 * t4
      t7 = t2 * t5 + t6
      t8 = za(i1, i3)
      t9 = za(i2, i3)
      t10 = zb(i3, i2)
      t11 = za(i3, i4)
      t12 = t8 * t4
      t13 = t9 * t10
      t14 = t11 * t2
      t15 = t14 + t12 + t13
      if ( real(t15) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t14 = cg * sqrt(t15 ** 2) + t12 + t13 + t14
      t15 = zb(i2, i1)
      t16 = za(i1, i4)
      t17 = zb(i4, i2)
      t18 = t3 * t15
      t19 = t16 * t1
      t20 = t17 * t5 + t18 + t19
      t21 = 0.1q1 / t14
      t19 = -2 * t12 * t20 * t21 + t18 + t19
      t22 = 2 * t8 * t20
      t23 = t22 * t2 * t21 - t17 * t3
      t24 = t10 * t3 - t16 * t2
      t25 = t8 * t19
      t17 = t4 * (t9 * (-t22 * t10 * t21 + t16 * t17) + t25)
      t22 = 2 * t9 * t20
      t26 = t22 * t2 * t21 + t1 * t3
      t22 = -t22 * t4 * t21 + t1 * t5
      t27 = 2 * t11 * t4 * t20 * t21 - t15 * t5
      t28 = t8 * (t22 * t10 + t19 * t4)
      t15 = 0.1q1 / t15
      t14 = 0.1q1 / t14
      t29 = 0.1q1 / t5
      t30 = 0.1q1 / t3
      t31 = 0.1q1 / t16
      t17 = 0.1q1 / t17
      t32 = 0.1q1 / t20
      t33 = t4 * t29
      t34 = t2 * t30
      t35 = t34 + t33
      t36 = t17 ** 2
      t37 = t17 * t36
      t38 = t9 ** 2
      t39 = t4 ** 2
      t40 = t4 * t39
      t41 = mt ** 2
      t42 = t1 * t15
      t11 = 0.1q1 / t11
      t27 = 0.1q1 / t27
      t28 = 0.1q1 / t28
      t43 = 0.1q1 / t23
      t44 = t1 * t24
      t45 = t2 * t19
      t46 = t4 * t23
      t47 = -t46 - t44 - t45
      t48 = t2 ** 2
      t49 = t9 * t24
      t50 = t49 * t29
      t33 = t33 * t9
      t51 = t45 * t43
      t52 = t39 * t9
      t53 = t4 * t15
      t54 = t12 * t2
      t55 = (t3 * t31 + t42) * t32
      t56 = t4 * t17
      t57 = t39 * t24
      t16 = t4 * (t39 * (t20 * t23 * t24 * t14 * t29 * t36 * t38 + t50 *
     # t25 * t42 * t20 * t14 * t36) + t4 * (-t50 * t42 * t20 * t14 * t17
     # + t19 * (t42 * t10 - t2) * t36 * t29 * t14 * t24 * t20 * t38) - t
     #55 * t2) + t29 * t15 * (t9 * (t4 * (-t56 * t24 * t26 + t11 * t7) -
     # t45 * t7 * (t11 * t43 + t56)) - t4 * t7 * (t8 * t2 * t28 + t27) *
     # t22 + t19 * (t57 * t36 * (t10 * t26 + t46) - t2 * (t10 * t5 + t16
     # * t4) * t11 ** 2 * t43) * t38) * t21 * t30 * t41
      t26 = t46 + t44
      t26 = -t23 * t26 + t50 * t26
      t58 = t50 - t23
      t59 = t43 ** 2
      t60 = t14 * t15 * t20
      t26 = t60 * t9 * (t29 * (t4 * (t44 * t43 + t4) + t34 * t49 * t43 *
     # (t43 * (t45 + t44) + t4) + t48 * t19 ** 2 * t59 + t44 * t45 * t59
     #) * t11 + t56 * t34 * (-t46 - t44 - t45) + t30 * t39 * (t12 * t26 
     #+ t13 * t26 + t45 * (t12 * t58 + t13 * t58)) * t36)
      t8 = t57 * t17 * t29 * t14 * t20 * t38 * (t46 * t19 * t36 * (t13 +
     # t12) + t60 * (-t46 * t38 * t10 ** 2 * t19 * t36 - t40 * t8 ** 2 *
     # t19 * t23 * t36 + t39 * t8 * t23 * t17 - t54 * t19 * t17 + t2 + t
     #13 * t17 * (t46 - t45)) * t30)
      t12 = -16
      t13 = 4
      t44 = 128
      ret = t44 * (t25 * t20 ** 2 * t14 ** 2 * t9 * t38 * t39 ** 2 * t10
     # * t24 * t23 * t30 * t29 * t15 * t37 + (t31 * (t29 * t6 + t2) + t4
     #2 * t35) * t32 * t21 * t7 * t41) + t13 * (t54 * t22 * t15 * t35 * 
     #t28 + t53 * t35 * t27 * t22 + t52 * (-t23 * t47 + t47 * t50) * t36
     # + t53 * (-t52 * t23 * t29 + t46 * (-t34 * t9 - t1) - t1 * t2 * (t
     #30 * t49 + t19)) * t17 + t15 * (t48 * t9 * t19 * t30 * t43 + t51 *
     # (t33 - t1) + t1 * t4 * (t43 * t50 - 1)) * t11) + 32 * t16 + t12 *
     # (-t18 * t38 * t40 * t19 * t24 * t23 * t29 * t37 + t5 * t48 * t32 
     #* (t30 * t42 + t31) + t3 * (t1 * t19 * t36 * t49 + t55) * t29 * t3
     #9 - t51 * t33 * t20 * t14 * t11 * t15) + 8 * t26 - 64 * t8

      hjetmass_qqbgg_triangle_pmpp_s124_mhsq_0 = ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     &  hjetmass_qqbgg_triangle_pmpp_s124_mhsq_0_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          real*16 cg

      t1 = za(i2, i3)
      t2 = zb(i3, i1)
      t3 = za(i1, i2)
      t4 = za(i2, i4)
      t5 = zb(i4, i3)
      t6 = t2 * t3 + t4 * t5
      t7 = za(i1, i3)
      t8 = zb(i3, i2)
      t9 = za(i3, i4)
      t10 = t7 * t2
      t11 = t1 * t8
      t12 = t9 * t5
      t13 = t11 + t12 + t10
      if ( real(t13) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t11 = cg * sqrt(t13 ** 2) + t10 + t11 + t12
      t12 = zb(i2, i1)
      t13 = za(i1, i4)
      t14 = zb(i4, i1)
      t15 = zb(i4, i2)
      t16 = t3 * t12
      t17 = t13 * t14
      t18 = t15 * t4 + t16 + t17
      t11 = 0.1q1 / t11
      t19 = 0.2q1 * t1
      t20 = -t19 * t2 * t18 * t11 + t14 * t4
      t21 = 0.2q1 * t9 * t2 * t18 * t11 - t12 * t4
      t10 = -0.2q1 * t10 * t18 * t11 + t16 + t17
      t16 = t7 * (t10 * t2 + t20 * t8)
      t17 = 0.2q1 * t7
      t22 = t2 * (t1 * (-t17 * t8 * t18 * t11 + t13 * t15) + t10 * t7)
      t15 = t17 * t18 * t5 * t11 - t15 * t3
      t17 = -t13 * t5 + t3 * t8
      t19 = t19 * t18 * t5 * t11 + t14 * t3
      t18 = 0.1q1 / t18
      t23 = 0.1q1 / t3
      t12 = 0.1q1 / t12
      t24 = 0.1q1 / t4
      t25 = 0.1q1 / t13
      t14 = t14 * t12
      t9 = 0.1q1 / t9
      t16 = 0.1q1 / t16
      t26 = 0.1q1 / t15
      t21 = 0.1q1 / t21
      t22 = 0.1q1 / t22
      t27 = t2 * t22
      ret = 0.128q3 * t11 * t18 * t6 * (t5 * (t14 * t23 + t25) + (t25 * 
     #t3 + t14) * t24 * t2) + 0.32q2 * t24 * t12 * t11 * t23 * (t1 * (t2
     # * (-t27 * t17 * t19 + t6 * t9) - t6 * (t26 * t9 + t27) * t10 * t5
     #) - t2 * t6 * (t5 * t7 * t16 + t21) * t20 + t10 * (-t5 * (t13 * t2
     # + t4 * t8) * t9 ** 2 * t26 + t2 ** 2 * t17 * t22 ** 2 * (t15 * t2
     # + t19 * t8)) * t1 ** 2)

      hjetmass_qqbgg_triangle_pmpp_s124_mhsq_0_rat = ret/32q0*(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_triangle_pmpp_s12_s123_0
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt
          real*16 cg

      t1 = za(i1, i2)
      t2 = za(i1, i3)
      t3 = zb(i3, i1)
      t4 = za(i2, i3)
      t5 = zb(i3, i2)
      t6 = t2 * t3
      t7 = t4 * t5
      t8 = t7 + t6
      if ( real(t8) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t8 = cg * sqrt(t8 ** 2) + t6 + t7
      t9 = za(i2, i4)
      t10 = za(i3, i4)
      t11 = zb(i2, i1)
      t12 = 0.1q1 / t8
      t13 = 2 * t1
      t14 = t13 * t10
      t15 = t11 * (t14 * t3 * t12 - t9)
      t16 = t1 * t11
      t17 = t16 * (-2 * t6 * t12 + 1)
      t18 = zb(i4, i1)
      t19 = zb(i4, i3)
      t20 = zb(i4, i2)
      t21 = t1 * (2 * t2 * t11 * t19 * t12 - t20)
      t22 = t3 * t17
      t13 = t7 * t6 * t13 * t11 * t12
      t23 = -t13 + t2 * (-t15 * t19 + t22)
      t24 = t10 * t21
      t13 = -t13 - t3 * (-t17 * t2 + t24)
      t25 = za(i1, i4)
      t26 = t25 * t18
      t27 = t9 * t20
      t14 = -t14 * t11 * t19 * t12 + t26 + t27
      t28 = 0.1q1 / t11
      t29 = 0.1q1 / t21
      t30 = 0.1q1 / t10
      t13 = 0.1q1 / t13
      t31 = 0.1q1 / t1
      t32 = 0.1q1 / t9
      t33 = t13 ** 2
      t34 = t13 * t33
      t35 = t29 ** 2
      t36 = t3 ** 2
      t37 = t3 * t36
      t38 = t5 ** 2
      t39 = t5 * t38
      t40 = t4 ** 2
      t41 = t18 ** 2
      t42 = t1 ** 2
      t43 = t7 * t13
      t44 = t42 * t40
      t45 = t17 * t13
      t46 = t28 * t31
      t47 = t28 * t32
      t48 = t4 * t19
      t49 = t13 * t21
      t50 = t28 * t18
      t51 = t50 * t25
      t52 = t20 * t28
      t53 = t19 * t13
      t54 = t17 * t19
      t55 = t1 * t18
      t56 = t55 * t5
      t57 = t47 * t54
      t58 = t49 * t4
      t59 = t27 + t26
      t60 = t21 ** 2
      t61 = t7 * t1
      t62 = t4 * t3
      t15 = 0.1q1 / t15
      t8 = 0.1q1 / t8
      t23 = 0.1q1 / t23
      t63 = t10 * t19
      t64 = -t63 - t6 - t7
      t65 = t63 + t6
      t66 = -t62 * t32 + t18
      t67 = t19 ** 2
      t68 = t54 * t35
      t69 = t3 * t29
      t70 = t3 * t21
      t71 = t1 * t32
      t72 = t53 * t3
      t73 = t8 * t17
      t74 = t1 * t8
      t24 = t4 * (t72 * t8 * (t22 * t58 + t55) * t5 + t37 * t19 * (t71 *
     # t12 * t23 + t73 * t21 * t33) * t2 + t73 * t67 * (t71 * t17 * t30 
     #* t35 + t24 * t36 * t33) + t74 * t36 * t4 * t21 * t33 * t66 * t38 
     #- t44 * t8 * t36 * t39 * t18 * t32 * t33) + t4 * (t8 * (t72 * (t70
     # + t54) + t71 * (t56 * (t69 + t68) + t36) * t30 + t36 * ((t63 + t6
     # + t7) * t3 * t60 + t5 * t1 * (t6 * t66 + t63 * t66) * t21 + t61 *
     # t54 * t32 * t64 - t42 * t4 * t38 * t18 * t32 * t65) * t33) + t36 
     #* (t15 * (t71 * t3 + t19) + t2 * t67 * t23) * t12)
      t55 = t32 * (t16 + t26) + t20
      t66 = t8 ** 2
      t72 = mt ** 2
      t73 = t22 * t16 * t21 * t32
      t75 = t16 * t32
      t22 = t22 * t21
      t2 = t36 * t4 * t1 * (-t73 * t66 * t4 * t40 * t39 * t34 + (t22 * t
     #55 * t34 * t8 + t75 * (t70 - t54) * t33 * t66) * t38 * t40 + t7 * 
     #(t66 * (t75 * (-t54 * t65 + t70 * t65) * t33 - t73 * (t67 * t10 **
     # 2 + t36 * t2 ** 2) * t34 + t75 * t53) + t22 * (t6 * t55 + t63 * t
     #55) * t34 * t8) - t72 * t3 * t12 ** 2 * t32 * (t2 * t19 * t23 + t1
     #5))
      t10 = -t51 - t1
      t15 = t13 * t4
      t22 = t5 * t1 * (2 * t4 * t11 * t19 * t12 + t18)
      t12 = t4 * (t47 * (t69 * (t22 + t54) * t30 + t13 * t36 * (t43 * t1
     #7 * (t14 * t19 - t22 - t70) + t54 + t22) - t68 * t7 * t14 * t30 **
     # 2) * t12 * t72 + t13 * t5 * t36 * (t71 * (t18 * (t74 * (t45 * t64
     # + 1) + t15 * (t54 * t8 + t70 * (-t45 - t8)) * t25) + t15 * t16 * 
     #t8 * (-t70 + t54)) + t15 * (t74 * t54 + t70 * (t45 * t10 - t74)) *
     # t20))
      t8 = t62 * t17 * (-t71 * t8 * t19 * t30 * t29 + t7 * t21 * (t32 * 
     #(t41 * t25 ** 2 * t28 + t11 * t42) + t9 * t20 ** 2 * t28) * t34 * 
     #t36 + t56 * (t10 * t32 - t52) * t33 * t3)
      t9 = -4
      ret = t9 * (t62 * (-t46 * t60 * t59 * t33 * t36 + t57 * t30 * t29 
     #+ (t54 * (-t46 * t21 * t59 + t7 * (t32 * (t51 + t1) + t52)) + t50 
     #* t5 * (t61 * t20 - t26 * t21)) * t33 * t3) + t36 * (t47 * t1 * t4
     #0 * t38 * t41 * t25 * t33 + t47 * t4 * t30 + t13 * (t44 * t38 * t1
     #3 * t32 + t21 * t28 - t43 * t21 * (t27 * t28 + t1)) * t18 + t49 * 
     #t48 * (t46 - t45)) + t50 * (t53 * (t7 + t17) + t30) * t3 + t57 * t
     #4 * t30 * (t56 + t54) * t35 + t54 * t28 * t30 * (-t48 * t31 + t18)
     # * t29 + t58 * (-t49 + t47 + t43 * (t32 * (t51 + t1) + t52)) * t37
     #) - 8 * t24 - 64 * t2 + 128 * t44 * t32 * t34 * t66 * t21 * t17 * 
     #t11 * t5 * t37 * (t63 * t7 + t6 * (t63 + t7)) - 32 * t12 + 16 * t8

      hjetmass_qqbgg_triangle_pmpp_s12_s123_0 = ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_triangle_pmpp_s12_s124_0
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt
          real*16 cg

      t1 = za(i1, i3)
      t2 = za(i2, i4)
      t3 = zb(i3, i1)
      t4 = zb(i4, i1)
      t5 = za(i1, i2)
      t6 = zb(i2, i1)
      t7 = za(i1, i4)
      t8 = zb(i4, i2)
      t9 = t7 * t4
      t10 = t2 * t8
      t11 = t9 + t10
      if ( real(t11) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t11 = cg * sqrt(t11 ** 2) + t10 + t9
      t12 = 0.1q1 / t11
      t13 = -2 * t9 * t12
      t14 = t5 * t6
      t15 = t14 * (1 + t13)
      t16 = zb(i3, i2)
      t17 = zb(i4, i3)
      t18 = t5 * (-2 * t7 * t6 * t17 * t12 - t16)
      t19 = za(i2, i3)
      t20 = za(i3, i4)
      t21 = t7 * t15
      t22 = t20 * t18
      t13 = t10 * t14 * t13
      t23 = t4 * (t21 + t22) + t13
      t24 = t6 * (-2 * t5 * t20 * t4 * t12 - t19)
      t13 = t7 * (t15 * t4 + t17 * t24) + t13
      t25 = 0.1q1 / t19
      t26 = 0.1q1 / t20
      t27 = 0.1q1 / t6
      t28 = 0.1q1 / t5
      t29 = 0.1q1 / t18
      t23 = 0.1q1 / t23
      t30 = t4 * t2
      t31 = t30 * t25
      t32 = t4 ** 2
      t33 = t18 ** 2
      t34 = t29 ** 2
      t35 = t23 ** 2
      t36 = t23 * t35
      t37 = t8 ** 2
      t38 = t17 ** 2
      t39 = t16 * t19
      t40 = t1 * t3
      t41 = t40 * t25
      t42 = t2 * t3 * t37
      t43 = t3 * t18
      t44 = t23 * t4
      t45 = t27 * (t3 + t31)
      t46 = t5 * t25
      t47 = t46 * t10
      t48 = t26 * t29
      t49 = t5 ** 2
      t50 = t15 * t23
      t51 = t8 * t3
      t52 = t15 * t17
      t53 = t4 * t18
      t24 = 0.1q1 / t24
      t13 = 0.1q1 / t13
      t11 = 0.1q1 / t11
      t54 = t2 ** 2
      t55 = t18 * t2
      t56 = t5 * t2
      t57 = t10 * t18
      t58 = t32 * t7
      t59 = t51 * t5
      t21 = t17 * t2 * (t4 * (-t58 * t46 * t12 * t13 + t4 * (-t42 * t49 
     #* t20 * t25 + t4 * t33 * t20 + (t43 * t20 + t31 * (-t22 + t21)) * 
     #t8 * t5) * t11 * t35 - t59 * t11 * t23) + t46 * t15 * t11 * (t59 -
     # t52) * t34 * t26) + t30 * (t11 * (t17 * (-t53 + t52) * t23 - t46 
     #* (t59 * t29 + t4) * t26 + t4 * (t15 * (t17 * (t46 * t54 * t37 - t
     #9 * t18 - t57) + t20 * (t47 - t18) * t38) + t8 * (t8 * (-t51 * t54
     # * t49 * t25 + t56 * ((-t5 * t3 * t7 - t55) * t25 * t4 + t43)) + t
     #53 * ((-t31 + t3) * t7 * t5 + t55)) + t58 * t33) * t35) + t4 * (t2
     #4 * (t46 * t4 - t17) + t7 * t38 * t13) * t12)
      t6 = t49 * t6 * t25
      t22 = t30 * t15 * (-t46 * t48 * t11 * t17 + t59 * (t27 * (-t41 - t
     #16) - t46) * t35 * t4 + t57 * (t27 * (t1 ** 2 * t3 ** 2 * t25 + t1
     #9 * t16 ** 2) + t6) * t36 * t32)
      t55 = t53 + t52
      t57 = mt ** 2
      t58 = t17 * t20
      t59 = t5 * t11
      t60 = t23 * t2
      t19 = t2 * (t27 * t25 * (t32 * (t44 * t18 - t26) + t52 * (t4 * (t4
     #4 - t48) + (-t32 * t18 * t35 + t29 * t26 ** 2) * t28 * (t1 * t4 + 
     #t19 * t8) * t2)) * t12 * t57 + t23 * t8 * t32 * (t46 * (t3 * (t59 
     #* (t50 * (t10 + t58 + t9) - 1) + t60 * (t52 * t11 + t53 * (t50 + t
     #11)) * t1) + t14 * t2 * t23 * t11 * t55) + t60 * (t59 * t52 + t53 
     #* (t59 + t50 * (t40 * t27 + t5))) * t16))
      t59 = t25 * (t14 + t40) + t16
      t60 = t11 ** 2
      t14 = t14 * t25
      t61 = t14 * t60
      t62 = t53 * t15
      t63 = t62 * t61
      t7 = t56 * t32 * (-t63 * t2 * t54 * t8 * t37 * t36 + t57 * t4 * t1
     #2 ** 2 * t25 * (t7 * t17 * t13 - t24) + (t62 * t59 * t11 * t36 + t
     #61 * t55 * t35) * t37 * t54 + t10 * (t36 * (t62 * (t58 * t59 + t9 
     #* t59) * t11 - t63 * (t38 * t20 ** 2 + t32 * t7 ** 2)) + t14 * (t5
     #8 * t55 + t9 * t55) * t60 * t35 - t61 * t17 * t23))
      t11 = -4
      ret = t11 * (t30 * (t27 * (-t51 * t17 * t23 + t52 * (-t40 * t10 * 
     #t4 * t35 + t48) * t25 + t44 * (t44 * t10 * t16 + (t50 * (t39 + t40
     #) - 1) * t28 * t17) * t18) + t53 * t8 * t35 * (t31 - t3) * t5 + t3
     #1 * t49 * t3 * t37 * t35) + t15 * (t17 * (-t44 * t3 * t27 + t2 * (
     #t18 - t10 * (t16 * t27 + t46)) * t35 * t32 + t48 * t3 * t27 * (t47
     # * t29 + 1)) + t48 * t2 * t27 * t28 * t38) + t4 * (t45 * t44 * t18
     # - t45 * t26 + t30 * (t43 * t27 * t8 * (t1 * (t31 - t3) - t39) + t
     #4 * (-1 - t27 * t28 * (t39 + t40)) * t33 + t42 * t27 * (t16 + t41)
     # * t5) * t35) - t2 * t15 ** 2 * t38 * t25 * t26 * t27 * t34) - 8 *
     # t21 + 16 * t22 + 32 * t19 - 64 * t7 + 128 * t6 * t60 * t36 * t18 
     #* t15 * t8 * t4 * t32 * t54 * (t58 * t9 + t10 * (t9 + t58))

      hjetmass_qqbgg_triangle_pmpp_s12_s124_0 = -ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     & hjetmass_qqbgg_triangle_pmpp_s34_mhsq_s12
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt
          real*16 cg

          cg = 1q0

      t1 = zb(i4, i3)
      t2 = za(i1, i3)
      t3 = zb(i3, i1)
      t4 = za(i2, i3)
      t5 = zb(i3, i2)
      t6 = za(i1, i4)
      t7 = zb(i4, i1)
      t8 = za(i2, i4)
      t9 = zb(i4, i2)
      t10 = t2 * t3
      t11 = t4 * t5
      t12 = t6 * t7
      t13 = t8 * t9
      t14 = t10 + t11 + t12 + t13
      t15 = za(i1, i2)
      t16 = za(i3, i4)
      t17 = zb(i2, i1)
      t14 = -4 * t15 * t16 * t17 * t1 + t14 ** 2
      t14 = cg * sqrt(t14) + t10 + t11 + t12 + t13
      t18 = (0.1q1 / 0.4q1)
      t19 = t14 ** 2
      t20 = t15 * t16
      t21 = t20 * t17 * t1
      t22 = t19 * t18 - t21
      t22 = 0.1q1 / t22
      t23 = t3 * t4 + t7 * t8
      t19 = t19 * t22
      t24 = t19 * t23
      t12 = t10 + t12
      t25 = t14 * t15
      t26 = t25 * t17 * t22
      t27 = (0.1q1 / 0.2q1)
      t28 = t19 * t18
      t29 = t28 * t15 * t17
      t30 = -t12 * t26 * t27 + t29
      t31 = t14 * t18
      t25 = t25 * t22
      t32 = t25 * (t27 * t2 * t17 * t1 - t31 * t9)
      t33 = t2 * t5 + t6 * t9
      t34 = t21 * t27 * t14 * t22
      t12 = t18 * t19 * t12 - t34
      t35 = -t25 * (t27 * t6 * t17 * t1 + t31 * t5)
      t36 = t27 * t4 * t17 * t1 + t31 * t7
      t37 = t14 * t16 * t22
      t38 = t30 * t12
      t39 = t35 * t37 * t36 + t38 - t24 * t26 * t33 / 8
      t40 = 2 * t15 * t17 + t14
      t21 = t27 * t14 ** 2 - 2 * t21
      t41 = 2 * t16 * t1 + t14
      t23 = t26 * t23
      t10 = t10 + t11
      t42 = t18 * t19 * t10 - t34
      t20 = t20 * t27
      t43 = t14 * t1 * t22
      t5 = t43 * (t20 * t5 + t31 * t6)
      t6 = t37 * t1
      t10 = t28 * t16 * t1 - t27 * t6 * t10
      t28 = t43 * (-t31 * t2 + t20 * t9)
      t37 = t20 * t7 + t31 * t4
      t22 = -t5 * t14 * t17 * t22 * t37 + t38 - t23 * t19 * t33 / 8
      t11 = t11 + t13
      t13 = -t27 * t26 * t11 + t29
      t11 = t18 * t19 * t11 - t34
      t18 = t25 * (-t27 * t8 * t17 * t1 + t31 * t3)
      t20 = t43 * (-t20 * t3 + t31 * t8)
      t26 = -t43 * t37
      t25 = t25 * t36
      t2 = t2 * t7 + t4 * t9
      t9 = t19 * t2
      t2 = t6 * t2
      t6 = 0.1q1 / t8
      t19 = 0.1q1 / t15
      t22 = 0.1q1 / t22
      t29 = 0.1q1 / t39
      t31 = 0.1q1 / t4
      t21 = 0.1q1 / t21
      t33 = 0.1q1 / t17
      t34 = t32 * t20
      t36 = t35 * t26
      t37 = t36 + t34
      t39 = t1 * t19
      t43 = t31 * t7
      t44 = t39 + t43
      t45 = -t20 * t7 + t26 * t3
      t46 = t30 + t42
      t47 = t10 + t11
      t48 = t41 ** 2
      t49 = t12 ** 2
      t50 = t12 * t49
      t51 = t1 ** 2
      t52 = t22 ** 2
      t53 = t22 * t52
      t54 = t23 ** 2
      t55 = t23 * t54
      t56 = t3 ** 2
      t57 = t21 ** 2
      t58 = t10 * t25
      t59 = t30 * t31
      t60 = t1 * t35
      t61 = t59 * t29
      t62 = t61 * t6
      t63 = t3 * t32
      t64 = t1 * t7
      t65 = t15 * t17
      t66 = t65 * t24
      t67 = t66 * t6
      t68 = t67 * t31
      t69 = t39 * t16 * t40
      t70 = t7 * t5
      t71 = t1 * t12
      t72 = t3 * t28
      t73 = t51 * t16
      t74 = t73 * t40
      t75 = t24 * t29
      t76 = t6 * t23
      t77 = -t12 - t10 - t11
      t78 = -t30 - t42 - t13
      t79 = t30 + t42 + t13
      t80 = t42 + t13
      t81 = t30 + t13
      t82 = t29 ** 2
      t83 = t29 * t82
      t84 = t24 ** 2
      t85 = t28 * t6
      t86 = t85 * t5
      t87 = t35 * t25
      t88 = t32 * t18
      t89 = t30 * t6
      t90 = t89 * t84
      t91 = t52 * t28 * t5
      t92 = t91 * t54
      t93 = t90 * t35
      t94 = t12 * t31
      t95 = t74 * t33 * t19
      t96 = t68 * t21 * t48 * (t92 * t77 + (-t34 * t10 + t36 * t77) * t8
     #2 * t30 * t24) + t31 * ((t33 * (t70 * t28 * t19 * t80 + (t70 * t78
     # + t72 * t81) * t6 * t12) + t86 * t79 * t21 * t41 * t24) * t52 * t
     #54 + t90 * t21 * t82 * t41 * (t12 * (t88 + t87) + t34 * t13)) * t4
     #0 * t16 * t1 + t95 * (-t93 * t82 * t80 + t94 * (-t30 - t13) * t52 
     #* t28 * t54)
      t97 = t12 + t10 + t11
      t98 = t3 * t12
      t99 = t1 * t28
      t100 = t16 * t1 * t40
      t101 = t100 * t6 * t31
      t102 = t101 * t21
      t103 = t21 * t31
      t104 = t41 * t11 * t21
      t105 = t100 * t21
      t106 = t29 * t6 * t84
      t107 = t106 * (-t88 * t100 * t57 * t41 * t31 + (t105 * t31 * (-t87
     # + t88 * (-1 + t104)) - t60 + (t65 * t37 * t31 + t60 * t97) * t21 
     #* t41) * t29 * t30)
      t108 = t88 + t87
      t109 = -t30 - t42 - t13
      t110 = t30 ** 2
      t111 = t28 ** 2
      t112 = t16 ** 2
      t113 = t40 ** 2
      t114 = t108 * t19
      t115 = t25 * t19
      t116 = t30 * t1
      t117 = t1 * t49
      t118 = t31 * t33
      t119 = (-t98 * t85 * t15 * t54 * t10 * t52 + t75 * (-t116 * t26 + 
     #t34 * t7)) * t21 * t41 - t3 * t54 * t111 * t52 - t90 * t65 * t34 *
     # t82 * (t12 + t11) * t57 * t48
      t120 = t114 * t51 * t84 * t110 * t6 * t82 * t31 * t57 * t33 * t113
     # * t112
      t90 = t105 * (t118 * t75 * (t7 * (t89 * t18 - t114) - t89 * t3 * t
     #25 + t116 * (t115 + t76)) + t90 * (-t60 * t30 * t19 * t33 + t103 *
     # (t88 * t10 + t36 * (t42 + t13)) * t41) * t82 + t118 * t52 * t54 *
     # (t19 * t111 * t3 * t109 - t117 * t109 * t6))
      t111 = t12 + t10 + t11
      t114 = t89 * t75
      t36 = (t28 * (t72 * t111 + t117) * t21 * t41 + t15 * t12 * t6 * (-
     #t70 + t72)) * t52 * t54 + t75 * t21 * (t36 * t41 * t7 + t114 * (-t
     #87 * t80 + t88 * (-t42 - t13)) * t33 * t21 * t19 * t113 * t112 * t
     #51)
      t34 = t31 * t36 + t31 * t119 + t54 * (t102 * t33 * (t70 - t72) * t
     #22 + t103 * (t41 * (t28 * (t71 * t47 - t70 * t97) + (-t1 * t50 + t
     #28 * (t17 * t24 * t5 - t98 * t11) + t70 * t12 * t47) * t6 * t15) +
     # t99 * (t70 * t30 * t19 + t98 * t42 * t6) * t33 * t40 * t16) * t52
     #) + t21 * t96 + t54 * (t74 * t28 * t19 * t31 * t21 * t33 * t22 + t
     #31 * (t28 * (-t71 * (t69 * t42 * t21 * t33 + 1) + t70) + t49 * (t1
     # + (-t1 * t47 + t70 - t72) * t21 * t41) * t6 * t15) * t52) + t75 *
     # (t33 * (t56 * t32 * t6 + t63 * t44 + t64 * (-t19 * t35 + t59)) + 
     #t41 * (t62 * t40 * t24 * t16 * t1 * (t35 * (t26 * t30 + t58) + t34
     # * t46) * t57 + t6 * (t59 * t15 * t45 - t60 * t24) * t21) + t68 * 
     #t37 * t57 * t48) - t76 * t71 * t15 * t41 * t21 * t24 * t31 * t22 +
     # t107 - t120 + t90
      t36 = 0.1q1 / t16
      t14 = 0.1q1 / t14
      t68 = t100 * t41 * t57
      t14 = mt ** 2 * t14
      t74 = t14 * t19
      t88 = t74 * t33
      t90 = t88 + t68
      t96 = t18 * t26 + t20 * t25
      t107 = t103 * t65
      t117 = t95 * t18
      t119 = t41 * t20
      t120 = t20 * t6
      t121 = t26 * t31
      t122 = t5 * t26
      t123 = t122 * t6
      t124 = t10 + t12
      t125 = t10 + t11
      t126 = t5 * t25
      t127 = t74 * t118
      t128 = t28 * t20
      t129 = t65 * t31
      t130 = t18 * t28
      t131 = t12 * t10
      t132 = t94 * t65
      t133 = t65 * t41
      t134 = t94 * t21
      t135 = t134 * t6
      t65 = -t88 * t94 * t86 * t52 + t134 * t86 * (t65 * t97 * t41 + t10
     #0 * t109) * t53
      t97 = t114 * t21 * (t107 * t48 * t20 * t26 + t117)
      t46 = t54 * (t101 * t57 * (t126 * t69 * t33 - t119 * t28) * t22 + 
     #t135 * (t133 * (t122 * (-t12 * t41 * t21 + 1) + t128) + t100 * ((t
     #28 * (t11 * t18 + t20 * (t30 + t13)) + t5 * (t26 * t46 + t58)) * t
     #21 * t41 - t130 - t126)) * t52) + t55 * t65 + t6 * (t54 * (t95 * t
     #21 * (t130 * t103 * t16 * t40 + t5) * t22 + t12 * (-t5 * (t127 * t
     #18 * t2 + t1) + t105 * (t39 * t33 * t5 * t78 + t103 * (t28 * (t124
     # * t18 + t20 * t42) + t126 * t12) * t41) + t129 * (-t122 * t125 + 
     #t128 * t77) * t57 * t48 + t1 * t5 * t125 * t21 * t41) * t52) + t55
     # * (-t91 * t71 * t16 * t40 * t41 * t31 * t57 + t132 * t57 * t28 * 
     #t5 * t48 * (t11 * (-t10 - t12) - t131) * t53) + t115 * t61 * t112 
     #* t51 * t113 * t57 * t24 * t18 * t33) - t103 * t72 * t69 * t23 * t
     #25 * t33 * t22 + t97
      t58 = t28 * t5
      t61 = t58 * t103 * t40 * t16
      t65 = t5 * t52 * t54
      t69 = t30 * t42
      t95 = t11 ** 2
      t97 = t10 ** 2
      t101 = t21 * t41
      t134 = t57 * t48
      t136 = t42 ** 2
      t137 = t13 ** 2
      t138 = t7 * t33
      t139 = t101 * t23
      t140 = t94 * t86 * t33 * t57 * t53 * t19 * t55 * t51 * (t137 + t11
     #0 + t136) * t113 * t112
      t43 = t101 * (t121 * t75 * t63 + (t1 * (-t24 * t5 * t6 - t94 * t26
     #) + t43 * (t122 + t128)) * t22 * t23)
      t2 = (t33 * (t72 * t44 + t85 * t56 + t64 * (-t19 * t5 + t94)) + t1
     #35 * t45 * t41 * t15) * t22 * t23 + t102 * (-t87 * t104 * t84 * t3
     #0 * t82 + t71 * t54 * t22 * t33 - t24 * t41 * t21 * (t130 + t126) 
     #* t22 * t23 + t91 * (t101 * t47 - 1) * t55) + t127 * t87 * t106 + 
     #t118 * t76 * t57 * t19 * t51 * (t75 * t108 + t92 * t109) * t113 * 
     #t112 + t132 * t86 * (1 + t134 * (t95 + t97)) * t53 * t55 + t43 - t
     #88 * t87 * t59 * t84 * t11 * t6 * t82 + t105 * (t75 * (-t115 * t11
     #8 * t63 + t76 * (-t103 * t37 * t41 + t60 * t19 * t33)) + t87 * t10
     #3 * t106 * t41 + t31 * t22 * t23 * (t18 * (-t138 * t28 * t19 + (t1
     #39 * t5 * t2 * t22 + t138) * t6 * t12) + t33 * (t19 * (t71 - t70) 
     #- t98 * t6) * t25)) + t140 + t134 * t129 * t76 * t22 * (t28 * (t20
     # * t24 + t65 * t50) + t122 * t24)
      t5 = t6 * t35
      t28 = t5 * t32
      t37 = t7 * t35
      t39 = t39 * t32
      t43 = t23 * t33
      t44 = t105 * t29
      t45 = t31 * t24
      t14 = t45 * (t139 * t99 * t22 + t5 * t33 * t19 * (t79 * t57 * t23 
     #* t113 * t32 * t112 * t51 + t14 * t30 * t20 * t9) * t82 * t24 + t9
     #3 * t21 * t32 * (t100 * (t101 * (t13 * t77 + t42 * t77 - t38) + t4
     #2 + t13) + t33 * t21 * t19 * (t13 * (t30 + t42) + t69) * t113 * t1
     #12 * t51 - t133 * t125) * t83) + t45 * (t6 * (t15 * (t101 * (-t38 
     #* t17 * t84 * t32 * t35 * t83 + t116 * t75 + t23 * t22 * (t70 - t7
     #2)) + t134 * t83 * t35 * t32 * t30 * t84 * t17 * (t11 * t124 + t13
     #1)) + t88 * t84 * t30 * t32 * t35 * t82) + t44 * (t84 * (t101 * t8
     #9 * t32 * t35 * t29 + t28 * t110 * (-t101 * t125 + 1) * t82) + t76
     # * t75 * t35 * t32 * (t101 * t77 + 1) + t43 * (t6 * (-t37 + t63) -
     # t39)))
      t38 = t1 * t110
      t45 = t32 * t19
      t17 = t75 * t17 * t32
      t47 = t32 ** 2
      t50 = t82 * t31 * t84
      t59 = t89 * t15
      t59 = t50 * (t3 * (t101 * t111 * t47 + t59 * t32) + t105 * (t89 * 
     #t63 * t42 + t37 * (-t89 * t13 + t45 * t81)) * t33 + t45 * t114 * t
     #35 * t57 * (t137 + t136) * t33 * t113 * t112 * t51 - t59 * t37)
      t5 = t59 + t50 * (-t3 * t47 + t134 * t67 * t35 * t32 * ((t97 + t49
     # + t95) * t29 * t30 - t12 - t10 - t11) + t105 * (t110 * (t1 * t80 
     #- t37 + t63) * t6 + t45 * (t63 * t78 - t38)) * t33 + t101 * (t32 *
     # (t116 * t10 + t37 * (-t10 - t11)) + t89 * (t37 * t111 + t63 * t77
     #) * t15)) + t31 * t29 * t84 * (t29 * (t32 * (-t116 + t37) + t89 * 
     #(t17 * t35 + t116) * t15) + t28 * t112 * t51 * t113 * t57 * t24 * 
     #t30 * t110 * t19 * t33 * t82 + t44 * (t33 * (t30 * (t6 * (t63 * t1
     #3 + t38) - t39 * t80) + t37 * t42 * (t45 - t89)) + t5 * t101 * (t3
     #0 * t20 * t9 + t24 * t32 * t80)) + t101 * (t32 * (t29 * (-t37 * t1
     #2 + t116 * (t12 + t11)) - t1) + (t35 * (t17 - t7) + t63 - t38 * t1
     #11 * t29) * t6 * t15))
      t9 = t18 * t6 + t25 * t31
      t15 = t94 * t76 * t96 * t22
      t17 = t24 * t9
      t3 = t33 * (-t3 * t7 * t36 + t74 * (t17 * t36 + t15)) + t41 * (t15
     # * t16 + t17) * t57 * t40 * t1 + t43 * t73 * t19 * t9 * t57 * t113
      t9 = t121 + t120
      t15 = -16
      t17 = 8
      ret = t27 * t5 + t15 * (t36 * (t33 * (t8 * t7 ** 2 * t31 + t4 * t5
     #6 * t6) + t66 * (t121 + t120) * t57 * t48) + t54 * (-t123 * t94 * 
     #t13 * t90 * t52 + t123 * t31 * t90 * t22) + t76 * t21 * t12 * (t11
     #5 * t103 * t51 * t112 * t113 * t18 * t33 + t117 + t119 * (t107 * t
     #41 * t26 - t1)) * t22 + t62 * t24 * (t68 * t96 + t88 * t96)) + t17
     # * (t21 * (t121 * t72 * t41 * t23 * t22 + t135 * t33 * t52 * t19 *
     # t54 * t113 * t112 * (t126 * t78 + t58 * (t13 * (-t30 - t42) - t69
     #) * t22 * t23 + t130 * t78) * t51 + t41 * (t12 * (t126 * t103 * t1
     #6 * t40 * t11 * t52 * t54 + t61 * (t10 * t79 + t11 * t79) * t53 * 
     #t55) + t49 * (t61 * t79 * t53 * t55 + t65) - t75 * t30 * t20) * t6
     # * t1) + t46) + 32 * t3 - 64 * t23 * (t40 * t1 * t57 * t41 * t9 + 
     #t88 * t9 * t36) + 2 * t34 - 4 * t2 + t14

      hjetmass_qqbgg_triangle_pmpp_s34_mhsq_s12 = ret/32q0/(0,1q0)
      return

      end function

      complex*32 function
     &  hjetmass_qqbgg_triangle_pmpp_s34_mhsq_s12_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          real*16 cg

          cg = 1q0

      t1 = za(i1, i3)
      t2 = zb(i3, i1)
      t3 = za(i2, i3)
      t4 = zb(i3, i2)
      t5 = za(i1, i4)
      t6 = zb(i4, i1)
      t7 = za(i2, i4)
      t8 = zb(i4, i2)
      t9 = t1 * t2
      t10 = t3 * t4
      t11 = t5 * t6
      t12 = t7 * t8
      t13 = t10 + t11 + t12 + t9
      t14 = za(i1, i2)
      t15 = za(i3, i4)
      t16 = zb(i2, i1)
      t17 = zb(i4, i3)
      t18 = t14 * t16
      t13 = -0.4q1 * t18 * t15 * t17 + t13 ** 2
      t13 = cg * sqrt(t13) + t10 + t11 + t12 + t9
      t19 = t3 * t2
      t20 = t7 * t6
      t21 = t20 + t19
      t22 = 0.1q1 / 0.2q1
      t23 = t22 * t13
      t24 = -t14 * (t5 * t16 * t17 + t23 * t4)
      t25 = t14 * (t3 * t16 * t17 + t23 * t6)
      t26 = 0.1q1 / 0.4q1
      t27 = t13 ** 2
      t28 = t18 * t15 * t17
      t29 = t26 * t27 - t28
      t30 = t11 + t9
      t31 = t15 * t2
      t32 = t30 * t16
      t33 = t27 * t14
      t19 = t20 + t19
      t20 = t5 ** 2 * t6 ** 2
      t34 = t1 ** 2 * t2 ** 2
      t35 = t4 * t19 * t1
      t19 = t8 * t19 * t5
      t36 = t18 * t13
      t37 = t17 ** 2
      t30 = t30 * t14
      t38 = t30 * t17
      t39 = t1 * t7
      t40 = t11 * t9 * t13
      t41 = t22 * t36 * ((t18 + t9 + t12) * t17 * t15 + t20 + t34 + t35 
     #+ t19)
      t42 = -t26 * t33 * (t31 * t8 + t32) + t18 * ((-t39 * t37 - t38) * 
     #t16 * t15 + t40) + t41
      t43 = t18 * t21
      t44 = t14 * t15
      t45 = t17 * (t23 * t5 + t44 * t4)
      t46 = t44 * t6
      t47 = -t17 * (t23 * t3 + t46)
      t48 = -t11 - t9
      t5 = t3 * t5
      t49 = t27 * t16
      t50 = t15 ** 2
      t32 = t32 * t17 * t15
      t4 = t4 * t6
      t19 = t22 * t36 * ((t18 + t10 + t11) * t17 * t15 + t20 + t34 + t35
     # + t19)
      t20 = t26 * t49 * (t14 * t48 + t5 * t17) + t18 * (t14 * (t4 * t17 
     #* t50 - t32) + t40) + t19
      t30 = -t26 * t49 * (t39 * t17 + t30) + t41 + t18 * (t14 * (-t2 * t
     #8 * t17 * t50 - t32) + t40)
      t4 = t26 * t33 * (t4 * t15 + t16 * t48) + t19 + t18 * ((t5 * t37 -
     # t38) * t16 * t15 + t40)
      t2 = t14 * (-t7 * t16 * t17 + t23 * t2)
      t5 = t11 + t9
      t9 = t22 * t13 * t5 - t28
      t19 = t23 * t18
      t5 = -t18 * t5 + t19
      t31 = t31 * t14
      t32 = t17 * (t23 * t7 - t31)
      t26 = t13 * t26
      t33 = t13 * t17
      t10 = t10 + t12
      t6 = t1 * t6 + t3 * t8
      t34 = t15 * t17
      t11 = t11 + t12
      t12 = 0.1q1 / t7
      t35 = 0.1q1 / t16
      t4 = 0.1q1 / t4
      t36 = 0.1q1 / t3
      t37 = 0.1q1 / t42
      t38 = 0.1q1 / t14
      t39 = t4 - t37
      t40 = 0.1q1 / t29 ** 2
      t41 = t2 * t47
      t40 = t38 * t36 * t35 * t40 * t12
      t42 = t40 * t5
      t48 = t37 ** 2
      t49 = t21 ** 2
      t50 = t4 ** 2
      t30 = 0.1q1 / t30
      t20 = 0.1q1 / t20
      t51 = -t30 + t20
      t52 = t20 ** 2
      t53 = t30 ** 2
      t18 = t47 * (-t18 * t10 + t19)
      t15 = 0.1q1 / t15
      t19 = -0.1q1 / t29
      t19 = t19 ** 2
      t10 = t25 * (t22 * t13 * t10 - t28)
      ret = -0.8q1 * t40 * t43 * t13 * (t9 * (t51 * t32 * t25 + t41 * t5
     #1) + (-t20 * t47 + t30 * t47 + t9 * (t25 * (-t34 * t11 + t34 * t23
     #) + t18) * t53 + t9 * (t2 * t34 * t6 - t18) * t52) * t45 * t43 + t
     #9 * t17 * (-t23 * t1 + t44 * t8) * (-t53 + t52) * t45 * t43 ** 2) 
     #- 0.64q2 * t43 * t38 * t15 * t35 * t19 * (t33 * (-t31 * t22 + t26 
     #* t7) * t12 - t33 * (t46 * t22 + t26 * t3) * t36) + 0.2q1 * t40 * 
     #t24 * t49 * t13 * t27 * (-t25 * t37 + t25 * t4 - t5 * (t47 * (t22 
     #* t13 * t11 - t28) + t10) * t48 + t10 * t5 * t50) + 0.16q2 * t21 *
     # t27 * t38 * t15 * t35 * t19 * (t12 * t2 + t25 * t36) + 0.4q1 * t4
     #2 * t21 * t27 * (t39 * t32 * t25 + t41 * t39) - t42 * t24 * t49 * 
     #t27 ** 2 * (-t6 * t32 * t50 + (-t50 + t48) * t14 * (t1 * t16 * t17
     # - t23 * t8) * t21)

      hjetmass_qqbgg_triangle_pmpp_s34_mhsq_s12_rat = ret/32q0*(0,1q0)
      return

      end function


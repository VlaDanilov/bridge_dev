      t1 = 0.1D1 / sman
      t2 = mH ** 2
      t3 = dilogc(-t1 * t2 + 1)
      t4 = dilogc(-t2 / tman + 1)
      t5 = 0.1D1 / uman
      t2 = dilogc(-t2 * t5 + 1)
      t6 = cdlogwrap(uman * t1)
      t5 = cdlogwrap(tman * t5)
      t1 = cdlogwrap(tman * t1)
      t7 = LogMums + LogMumt + LogMumu
      t8 = sman * tman
      t9 = sman + tman
      t10 = sman ** 2
      t11 = sman * t10
      t12 = tman ** 2
      t13 = tman * t12
      t14 = t8 * t9
      t9 = (t9 * uman + t10 + t12) * uman + t14
      t15 = sman + tman
      t16 = uman * t15
      t15 = t8 * t15
      t17 = (t16 + t10 + t12) * uman + t15
      t18 = uman ** 3
      t19 = t13 + t18 + t11
      t1 = t1 ** 2
      t5 = t5 ** 2
      t6 = t6 ** 2
      t20 = t8 * uman
      t21 = t13 + t18 + t11
      t21 = LogMums * t21 + LogMumt * t21 + LogMumu * t21
      t22 = pi ** 2
      t23 = t22 - t1 - t5 - t6
      t24 = -t3 - t4 - t2
      t25 = t19 * LogMuMtop
      t26 = t20 * cA
      t27 = t9 * LogMumt
      t28 = t9 * LogMumu
      t29 = LogMums * t9
      t30 = (0.7D1 / 0.30D2)
      t31 = (0.7D1 / 0.60D2)
      t1 = t20 * (-cA * t22 / 12 - (0.2D1 / 0.5D1) * LogMuMtop * tr) + (
     #0.289D3 / 0.1080D4) * t19 * cA + (0.881D3 / 0.1080D4) * cA * t9 + 
     #(0.77D2 / 0.360D3) * cA * (t29 + t27 + t28) + (0.77D2 / 0.1080D4) 
     #* cA * t21 + (0.61D2 / 0.45D2) * cF * t9 + (0.61D2 / 0.135D3) * t1
     #9 * cF + (0.7D1 / 0.120D3) * cA * (t1 * t9 + t22 * ((uman * (-sman
     # - tman) - t10 - t12) * uman - t14) + t5 * t9 + t6 * t9) + (0.7D1 
     #/ 0.180D3) * t25 + t31 * (Kg * t19 + t17 * LogMuMtop + t20 * (t1 +
     # t5 + t6) * cA) + (0.7D1 / 0.20D2) * Kg * t9 + t30 * (cA * (t17 * 
     #t2 + t17 * t3 + t17 * t4) + ((-t16 - t10 - t12) * uman - t15) * Lo
     #gMuMtop * tr) - (0.2D1 / 0.3D1) * t8 * tr * uman * t7 - (0.7D1 / 0
     #.360D3) * cA * (t11 * t23 + t13 * t23 + t18 * t23) - (0.7D1 / 0.18
     #D2) * tr * (t29 + t18 + t11 + t13 + t27 + t28) - (0.7D1 / 0.54D2) 
     #* tr * t21 - (0.35D2 / 0.27D2) * tr * t9
      t5 = -(0.7D1 / 0.90D2)
      c1mt2 = t5 * (cA * (t11 * t24 + t13 * t24 + t18 * t24) + t25 * tr)
     # + (0.7D1 / 0.15D2) * t26 * (t3 + t4 + t2) + (0.13D2 / 0.30D2) * t
     #26 * t7 - t20 * LogMuMtop * (cA - 1) / 5 + t1 + t20 * ((0.3D1 / 0.
     #5D1) * Kg + (0.173D3 / 0.100D3) * cA + (0.37D2 / 0.15D2) * cF - (0
     #.44D2 / 0.15D2) * tr)

